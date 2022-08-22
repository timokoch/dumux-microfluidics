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
 * \brief Compute volume of intersections
 */
#ifndef DUMUX_MICROFLUIDIC_COMPUTE_VOLUME_HH
#define DUMUX_MICROFLUIDIC_COMPUTE_VOLUME_HH

#include <array>
#include <vector>
#include <cmath>
#include <utility>

#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>

namespace Dumux {

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
double computeReservoirVolume(const std::vector<std::array<Point, 4>>& intersections)
{
    double volume = 0.0;
    for (const auto& tet : intersections)
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
auto computeReservoirVolumeWithSplit(const std::vector<std::array<Point, 4>>& intersections, const double splitY)
{
    double vol0 = 0.0, vol1 = 0.0;
    constexpr std::size_t numCorners = 4;
    for (const auto& tet : intersections)
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

} // end namespace Dumux

#endif
