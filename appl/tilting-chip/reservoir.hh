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

namespace Dumux::Microfluidic {

class Reservoir
{
    using ALU = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;
    using GridView = ALU::LeafGridView;
    using EntitySet = Dumux::GridViewGeometricEntitySet<GridView, 0>;
    using Point = Dune::FieldVector<double, 3>;

public:
    // construct a grid representation of the reservoir geometry as
    // a triangulation (volume represented by a collections of tetrahedrons)
    // This offers flexibility concerning the geometry. However, for computing
    // intersections it might be more efficient to realize this via a CAD kernel
    // representation, e.g. in OpenCascade
    Reservoir()
    {
        gridManager_.init();
        const auto& gridView = gridManager_.grid().leafGridView();

        // write reservoir geometry to VTK file for debugging
        Dune::VTKWriter<GridView> vtkWriter(gridView);
        vtkWriter.write("grid", Dune::VTK::base64);

        // compute axis-aligned bounding box tree representation
        aabbTree_ = std::make_shared<Dumux::BoundingBoxTree<EntitySet>>(
            std::make_shared<EntitySet>(gridView)
        );

        // compute volume and bounding box
        for (const auto& element : elements(gridView))
        {
            const auto geometry = element.geometry();
            reservoirVolume_ += geometry.volume();

            for (int i = 0; i < geometry.corners(); ++i)
            {
                const auto& corner = geometry.corner(i);
                for (int dir = 0; dir < 3; ++dir)
                {
                    lowerLeft_[dir] = std::min(lowerLeft_[dir], corner[dir]);
                    upperRight_[dir] = std::max(upperRight_[dir], corner[dir]);
                }
            }
        }
    }

    // axis-aligned bounding box tree representation of the reservoir geometry
    const Dumux::BoundingBoxTree<EntitySet>& boundingBoxTree() const
    { return *aabbTree_; }

    // the bounding box enclosing all of the geometry
    std::pair<Point, Point> boundingBox() const
    { return std::make_pair(lowerLeft_, upperRight_); }

    // the reservoir volume
    double volume()
    { return reservoirVolume_; }

private:
    Dumux::GridManager<ALU> gridManager_;
    std::shared_ptr<Dumux::BoundingBoxTree<EntitySet>> aabbTree_;
    Point lowerLeft_{1e100}, upperRight_{-1e100};
    double reservoirVolume_{0.0};
};

} // end namespace Dumux::Microfluidic

#endif
