// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Reservoir geometry representation
 */
#ifndef DUMUX_MICROFLUIDIC_RESERVOIR_HH
#define DUMUX_MICROFLUIDIC_RESERVOIR_HH

#include <array>
#include <vector>
#include <cmath>
#include <utility>

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/nonlinear/findscalarroot.hh>

#include "intersections.hh"

namespace Dumux::Microfluidic {

/*!
 * \brief Compute the volume of a tetrahedron
 */
template<class Point>
double tetrahedronVolume(const std::array<Point, 4>& corners)
{
    using std::abs;
    return 1.0/6.0 * abs(
        Dumux::tripleProduct(
            corners[3]-corners[0],
            corners[1]-corners[0],
            corners[2]-corners[0]
        )
    );
}

/*!
 * \brief Compute the total volume of an intersection triangulation
 */
template<class Point>
double computeReservoirVolume(const std::vector<std::array<Point, 4>>& triangulation)
{
    double volume = 0.0;
    for (const auto& tet : triangulation)
    {
        const auto vol = tetrahedronVolume(tet);
        volume += std::isfinite(vol) ? vol : 0; // maybe there are broken intersections
    }
    return volume;
}

/*!
 * \brief Compute the total volume of an intersection triangulation but split into two
 * splitY is the split plane position and is assumed to be orthogonal to the y-axis
 * we assume that the reservoir is oriented such that a split in y corresponds to a
 * a split through the middle of the reservoir and each part belongs to one channel
 */
template<class Point>
auto computeReservoirVolumeWithSplit(const std::vector<std::array<Point, 4>>& triangulation, const double splitY)
{
    double vol0 = 0.0, vol1 = 0.0;
    constexpr std::size_t numCorners = 4;
    for (const auto& tet : triangulation)
    {
        const auto vol = tetrahedronVolume(tet);
        if (std::isfinite(vol)) // maybe there are broken intersections
        {
            for (int i = 0; i < numCorners; ++i)
            {
                if (tet[i][1] < splitY)
                    vol0 += vol/numCorners;
                else
                    vol1 += vol/numCorners;
            }
        }
    }
    return std::pair{ vol0, vol1 };
}

/*!
 * \brief Reservoir representation
 */
class Reservoir
{
    using Point = Dune::FieldVector<double, 3>;

    // internal grid representation of the reservoir geometry
    using Grid = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;
    using GridView = Grid::LeafGridView;
    using EntitySet = Dumux::GridViewGeometricEntitySet<GridView, 0>;

    // geometry type for helper hexahedron for the water table calculation
    template <class ct>
    struct HexTraits : public Dune::MultiLinearGeometryTraits<ct>
    {
        template<int mydim, int cdim>
        struct CornerStorage { using Type = std::array<Point, 8>; };

        template<int mydim>
        struct hasSingleGeometryType
        {
            static const bool v = true;
            static const unsigned int topologyId = Dune::GeometryTypes::cube(mydim).id();
        };
    };

    using HexGeometry = Dune::MultiLinearGeometry<double, 3, 3, HexTraits<double>>;

public:
    // default constructor (geometry is read from input file "Grid.File")
    Reservoir() { initialize_(); }

    // the reservoir volume (i.e. total volume of the geometry)
    double volume() const
    { return reservoirVolume_; }

    // information about a fluid body inside the reservoir
    struct FluidBody
    {
        // a triangulation of the fluid body (collection of tetrahedra)
        std::vector<std::array<Point, 4>> triangulation;

        // the total volume of the fluid body
        double volume;

        // the water table (plane representation in terms of a point and a normal vector)
        Point waterTablePoint;
        Point waterTableNormal;
    };

    // compute representation of the fluid body at a given set of rotation angles (pitch, roll) and fluid volume
    FluidBody computeFluidBody(const std::array<double, 2>& angle, double fluidVolume) const
    {
        // create helper hexahedron
        auto hexahedronCorners = hexahedronCorners_;

        // adjust the top plane of the helper hexahedron to simulate the water table angle
        Point rotationCenter(0.0);
        for (int i = 4; i < 8; ++i)
            rotationCenter += hexahedronCorners[i];
        rotationCenter /= 4.0;

        for (int i = 4; i < 8; ++i)
        {
            const auto heightDiffX = std::tan(angle[0])*(hexahedronCorners[i][1] - rotationCenter[1]);
            const auto heightDiffY = std::tan(angle[1])*(hexahedronCorners[i][0] - rotationCenter[0]);
            hexahedronCorners[i][2] += heightDiffX + heightDiffY;
        }

        // create a geometry (needed for the intersection interface)
        auto hexCut = HexGeometry(Dune::GeometryTypes::cube(3), hexahedronCorners);

        // we store the results in this object
        FluidBody fluidBody;

        // compute intersections of the helper hexahedron with the reservoir geometry
        auto treeIntersections = Dumux::intersectingEntities(hexCut, boundingBoxTree_);
        fluidBody.triangulation = Dumux::convertIntersections<Point>(treeIntersections);

        // we regularize the fluid volume to make sure the bracket root finding algorithm works
        fluidVolume = std::clamp(fluidVolume, 1e-2, reservoirVolume_-1e-2);

        // find the water table height by optimizing the reservoirVolume to the given one
        // formally we minimize the difference between the given volume and the intersection
        // volume of the intersection body between the reservoir and the half-space
        // below the water table plane
        auto localCorners = hexahedronCorners; // avoid reallocation every iteration
        const auto residual = [&](const double h)
        {
            // shift water table height
            localCorners = hexahedronCorners;
            for (int i = 4; i < 8; ++i)
                localCorners[i][2] += h;

            hexCut = HexGeometry(Dune::GeometryTypes::cube(3), localCorners);
            Dumux::intersectingEntities(treeIntersections, hexCut, boundingBoxTree_);
            Dumux::convertIntersections<Point>(fluidBody.triangulation, treeIntersections);
            const auto iVol = computeReservoirVolume<Point>(fluidBody.triangulation);
            return iVol - fluidVolume;
        };

        // water table height is the local variable
        // Brent's method is due to Brent 1971: https://doi.org/10.1093/comjnl/14.4.422
        double h = 0.0;
        try { h = Dumux::findScalarRootBrent(-11.0, 5.0, residual, 1e-4, 2000); }
        catch (const Dune::InvalidStateException& e)
        { h = Dumux::findScalarRootBrent(-11.1, 5.1, residual, 1e-4, 2000); }

        // store the resulting volume and intersections
        // (each application of "residual" stores the intersections as side effect)
        fluidBody.volume = residual(h) + fluidVolume;

        localCorners = hexahedronCorners;
        for (int i = 4; i < 8; ++i)
            localCorners[i][2] += h;

        const auto ab = localCorners[5] - localCorners[4];
        const auto ac = localCorners[6] - localCorners[4];

        // store the water table plane information
        fluidBody.waterTableNormal = Dumux::crossProduct(ab, ac);
        fluidBody.waterTableNormal /= fluidBody.waterTableNormal.two_norm();
        fluidBody.waterTablePoint = localCorners[4];

        return fluidBody;
    }

private:
    Dumux::GridManager<Grid> gridManager_;
    Dumux::BoundingBoxTree<EntitySet> boundingBoxTree_;
    double reservoirVolume_{0.0};

    // helper construct for intersecting the reservoir geometry with a hexahedron
    // the top face of the hexahedron is aligned with the water table
    // and the intersection body is the fluid-filled part of the reservoir
    std::array<Point, 8> hexahedronCorners_;

    // initialize internal geometry representation and helpers
    // for the computation of the water table
    //
    // we construct a grid representation of the reservoir geometry as
    // a triangulation (volume represented by a collections of tetrahedrons)
    // This offers flexibility concerning the geometry. However, for computing
    // intersections it might be more efficient to realize this via a CAD kernel
    // representation, e.g. in OpenCascade
    void initialize_()
    {
        gridManager_.init();
        const auto& gridView = gridManager_.grid().leafGridView();

        // write reservoir geometry to VTK file for debugging
        Dune::VTKWriter<GridView> vtkWriter(gridView);
        vtkWriter.write("grid", Dune::VTK::base64);

        // compute axis-aligned bounding box tree representation
        boundingBoxTree_ = Dumux::BoundingBoxTree<EntitySet>(
            std::make_shared<EntitySet>(gridView)
        );

        // compute volume and bounding box
        Point lowerLeft(1e100), upperRight(-1e100);
        for (const auto& element : elements(gridView))
        {
            const auto geometry = element.geometry();
            reservoirVolume_ += geometry.volume();

            for (int i = 0; i < geometry.corners(); ++i)
            {
                const auto& corner = geometry.corner(i);
                for (int dir = 0; dir < 3; ++dir)
                {
                    lowerLeft[dir] = std::min(lowerLeft[dir], corner[dir]);
                    upperRight[dir] = std::max(upperRight[dir], corner[dir]);
                }
            }
        }

        // create helper hexahedron for the water table calculation
        const auto height = upperRight[2]-lowerLeft[2];
        lowerLeft[2] -= height;
        upperRight[2] += height;

        hexahedronCorners_ = std::array<Point, 8>{
            lowerLeft,
            { upperRight[0], lowerLeft[1], lowerLeft[2] },
            { lowerLeft[0], upperRight[1], lowerLeft[2] },
            { upperRight[0], upperRight[1], lowerLeft[2] },
            { lowerLeft[0], lowerLeft[1], upperRight[2] },
            { upperRight[0], lowerLeft[1], upperRight[2] },
            { lowerLeft[0], upperRight[1], upperRight[2] },
            upperRight,
        };

        for (int i = 4; i < 8; ++i)
            hexahedronCorners_[i][2] -= 1.0*height;
    }
};

} // end namespace Dumux::Microfluidic

#endif
