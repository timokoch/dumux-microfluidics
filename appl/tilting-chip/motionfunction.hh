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
 * \brief The motion function of the tilting platform
 */
#ifndef DUMUX_MICROFLUIDIC_MOTION_FUNCTION_HH
#define DUMUX_MICROFLUIDIC_MOTION_FUNCTION_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

namespace Dumux::Microfluidic {

/*!
 * \brief Base class for tilting function
 */
class MotionFunction
{
public:
    virtual ~MotionFunction() = default;
    // rotation angles at time t
    virtual std::pair<double, double> rotationAngles(double t) const = 0;
};

/*!
 * \brief Tilting function that models rotation with constant tilt
 */
class ConstantTiltMotionFunction : public MotionFunction
{
public:
    explicit ConstantTiltMotionFunction(double rotationsPerSecond)
    : rotationsPerSeconds_(rotationsPerSecond)
    {
        // we use a constant tilt angle for this motion function
        const auto angleDegree = getParam<double>("Problem.Angle", 18.0);
        tiltAngle_ = angleDegree/180.0*M_PI;
        sinTiltAngle_ = std::sin(tiltAngle_);
    }

    // rotation angles (pitch, roll) at time curTime
    std::pair<double, double> rotationAngles(double curTime) const final
    {
        const auto t = curTime*rotationsPerSeconds_*2.0*M_PI;
        const auto sinT = std::sin(t);
        const auto cosT = std::cos(t);
        const auto gamma = std::asin(-sinTiltAngle_*sinT);
        const auto beta = std::asin(-sinTiltAngle_*cosT/std::sqrt(1.0 - sinTiltAngle_*sinTiltAngle_*sinT*sinT));
        return std::make_pair(beta, gamma);
    }

private:
    double tiltAngle_, sinTiltAngle_;
    double rotationsPerSeconds_;
};

} // end namespace Dumux::Microfluidic

#endif
