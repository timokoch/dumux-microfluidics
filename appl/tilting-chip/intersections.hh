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
 * \brief Algorithms that finds which geometric entities intersect
 */
#ifndef DUMUX_MICROFLUIDIC_GEOMETRY_INTERSECTIONS_HH
#define DUMUX_MICROFLUIDIC_GEOMETRY_INTERSECTIONS_HH

#include <array>
#include <vector>

#include <test/geometry/writetriangulation.hh>
#include <dumux/geometry/intersectingentities.hh>

namespace Dumux {

/*!
 * \brief Computes all entities that intersect with the given geometry
 *
 * In contrast to the implementation in Dumux we provide the container for the intersections
 * so we don't have to reallocate memory every time
 */
template<class Geometry, class EntitySet>
inline void
intersectingEntities(std::vector<IntersectionInfo<Geometry::coorddimension, typename Geometry::ctype, typename EntitySet::ctype>>& intersections,
                     const Geometry& geometry,
                     const BoundingBoxTree<EntitySet>& tree)
{
    // check if the world dimensions match
    static_assert(int(Geometry::coorddimension) == int(EntitySet::dimensionworld),
        "Can only intersect geometry and bounding box tree of same world dimension");

    // Create data structure for return type
    intersections.clear();
    using ctype = typename IntersectionInfo<Geometry::coorddimension, typename Geometry::ctype, typename EntitySet::ctype>::ctype;
    static constexpr int dimworld = Geometry::coorddimension;

    // compute the bounding box of the given geometry
    std::array<ctype, 2*Geometry::coorddimension> bBox;
    ctype* xMin = bBox.data(); ctype* xMax = xMin + Geometry::coorddimension;

    // Get coordinates of first vertex
    auto corner = geometry.corner(0);
    for (std::size_t dimIdx = 0; dimIdx < dimworld; ++dimIdx)
        xMin[dimIdx] = xMax[dimIdx] = corner[dimIdx];

    // Compute the min and max over the remaining vertices
    for (std::size_t cornerIdx = 1; cornerIdx < geometry.corners(); ++cornerIdx)
    {
        corner = geometry.corner(cornerIdx);
        for (std::size_t dimIdx = 0; dimIdx < dimworld; ++dimIdx)
        {
            using std::max;
            using std::min;
            xMin[dimIdx] = min(xMin[dimIdx], corner[dimIdx]);
            xMax[dimIdx] = max(xMax[dimIdx], corner[dimIdx]);
        }
    }

    // Call the recursive find function to find candidates
    intersectingEntities(geometry, tree,
                         bBox, tree.numBoundingBoxes() - 1,
                         intersections);
}

namespace Detail {

template <std::size_t N, class T, std::size_t... I>
constexpr std::array<T, N> vecToArrayImpl(const std::vector<T>& vec, std::index_sequence<I...>)
{ return { {vec[I]...} }; }

template <std::size_t N, class T>
constexpr std::array<T, N> vecToArray(const std::vector<T>& vec)
{ assert(vec.size() == N); return vecToArrayImpl<N>(vec, std::make_index_sequence<N>{}); }

} // end namespace Detail

/*!
 * \brief Convert intersection format to simple corner storage
 * This returns a vector of simplices (given by a vector of points)
 */
template<class Point, class TreeIntersections>
void convertIntersections(std::vector<std::array<Point, 4>>& intersections, const TreeIntersections& treeIntersections)
{
    intersections.clear();
    intersections.reserve(treeIntersections.size());
    for (const auto& is : treeIntersections)
    {
        assert(is.corners().size() == 4);
        intersections.emplace_back(Detail::vecToArray<4>(is.corners()));
    }
}

/*!
 * \brief Convert intersection format to simple corner storage
 * This returns a vector of simplices (given by a vector of points)
 */
template<class Point, class TreeIntersections>
auto convertIntersections(const TreeIntersections& treeIntersections)
{
    std::vector<std::array<Point, 4>> intersections;
    convertIntersections(intersections, treeIntersections);
    return intersections;
}

/*!
 * \brief Write intersections to file for debugging and visualization
 */
template<class Point>
void writeIntersections(const std::vector<std::array<Point, 4>>& intersections, const std::string& name)
{ Dumux::writeVTUTetrahedron(intersections, name); }

} // end namespace Dumux

#endif
