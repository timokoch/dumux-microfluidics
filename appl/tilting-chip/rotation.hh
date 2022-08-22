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
 * \brief Create rotation functions
 */
#ifndef DUMUX_MICROFLUIDIC_ROTATION_HH
#define DUMUX_MICROFLUIDIC_ROTATION_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>

namespace Dumux {

/*!
 * \brief Create a function that rotated a point around an axis with a given angle
 */
template<class ctype>
auto make3DRotation(const Dune::FieldVector<ctype, 3>& rotationAxis,
                    const ctype rotationAngle)
{
    using std::sin; using std::cos;
    const ctype sinAngle = sin(rotationAngle);
    const ctype cosAngle = cos(rotationAngle);
    return [=](Dune::FieldVector<ctype, 3> p){
        auto tp = p;
        tp *= cosAngle;
        tp.axpy(sinAngle, Dumux::crossProduct({rotationAxis}, p));
        tp.axpy((1.0-cosAngle)*(rotationAxis*p), rotationAxis);
        return tp;
    };
}

} // end namespace Dumux

#endif
