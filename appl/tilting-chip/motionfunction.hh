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
#include <vector>
#include <algorithm>

#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/io/container.hh>

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
    std::pair<double, double> rotationAngles(double curTime) const override
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

/*!
 * \brief Tilting function that models rotation with constant tilt but variable speed
 */
class LinearInterpolationSpeedConstantTiltMotionFunction : public ConstantTiltMotionFunction
{
public:
    explicit LinearInterpolationSpeedConstantTiltMotionFunction(double rotationsPerSecond)
    : ConstantTiltMotionFunction(rotationsPerSecond)
    , rotationsPerSeconds_(rotationsPerSecond)
    {
        x_ = getParam<std::vector<double>>("Problem.MotionPointsX");
        y_ = getParam<std::vector<double>>("Problem.MotionPointsY");
        if (x_.size() != y_.size())
            DUNE_THROW(ParameterException, "MotionPointsX size != MotionPointsY size");
        std::sort(x_.begin(), x_.end());
        std::sort(y_.begin(), y_.end());
        if (x_[0] > 1e-20 || y_[0] > 1e-20)
            DUNE_THROW(ParameterException, "First control point has to be (0, 0)");
        if (x_.back() < 0.5 || y_.back() < 0.5)
            DUNE_THROW(ParameterException, "Last control point has to be (0.5, 0.5)");

        const auto plotTime = Dumux::linspace(0.0, 1.0/rotationsPerSeconds_, 100);
        std::vector<Dune::FieldVector<double, 5>> output; output.reserve(100);
        for (int i = 0; i < 100; ++i)
        {
            const auto [a, b] = rotationAngles(plotTime[i]);
            const auto t = mapToUnitInterval_(plotTime[i]);
            output.push_back(Dune::FieldVector<double, 5>{{ t, transform_(t), plotTime[i], a, b }});
        }

        writeContainerToFile(output, "angles.txt", 10);
    }

    // rotation angles (pitch, roll) at time curTime
    std::pair<double, double> rotationAngles(double curTime) const override
    {
        // get value in [0, 1)
        const auto t = mapToUnitInterval_(curTime);

        // apply the a transformation of the parameter-space
        // (takes a number between 0 and 0.5 and spits out a number between 0 and 0.5)
        // tt is in [0, 1)
        const auto tt = transform_(t);

        // convert to time and forward to the usual implementation
        return ConstantTiltMotionFunction::rotationAngles(tt/rotationsPerSeconds_);
    }

private:
    double mapToUnitInterval_(double curTime) const
    {
        return curTime*rotationsPerSeconds_ - std::floor(curTime*rotationsPerSeconds_);
    }

    double transform_(const double tOne) const
    {
        // we take a number between 0 and 0.5 because we want the
        // sequence to be the same for half of time (left, right channel)
        const auto offset = tOne > 0.5 ? 0.5 : 0.0;
        const auto t = tOne - offset;

        if (t <= x_[0])
            return y_[0];
        if (t >= x_.back())
            return y_.back();

        const auto lookUpIndex = std::distance(x_.begin(), std::lower_bound(x_.begin(), x_.end(), t));
        if (lookUpIndex < 1)
            DUNE_THROW(Dune::InvalidStateException, "Wrong interpolation in LinearInterpolationSpeedConstantTiltMotionFunction");
        const auto gamma = (t - x_[lookUpIndex-1])/(x_[lookUpIndex] - x_[lookUpIndex-1]); // gamma in [0, 1)
        const auto tt = y_[lookUpIndex-1] + gamma*(y_[lookUpIndex] - y_[lookUpIndex-1]);

        // go back to [0, 1)
        return tt + offset;
    }

    double rotationsPerSeconds_;
    std::vector<double> x_, y_;
};

} // end namespace Dumux::Microfluidic

#endif
