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
 * \brief The channel state
 */
#ifndef DUMUX_MICROFLUIDIC_CHANNEL_HH
#define DUMUX_MICROFLUIDIC_CHANNEL_HH

#include <iostream>
#include <string>
#include <array>

#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/io/format.hh>

#include "rotation.hh"
#include "reservoir.hh"

namespace Dumux::Microfluidic {

class ChannelState
{
public:
    ChannelState(double pressure, double availableFluidVolume, bool isDry)
    : pressure_(pressure)
    , availableFluidVolume_(availableFluidVolume)
    , isDry_(isDry)
    {}

    double pressure() const
    { return pressure_; }

    double availableFluidVolume() const
    { return availableFluidVolume_; }

    bool isDry() const
    { return isDry_; }

private:
    double pressure_;
    double availableFluidVolume_;
    bool isDry_;
};

// compute channel states of a reservoir
// we assume that there are two channels
// the channel inlet/outlet position is given by the parameters
// "Problem.MeasurementPoint1"/"Problem.MeasurementPoint2" in the input file
auto computeChannelStates(const Microfluidic::Reservoir& reservoir,
                          const std::array<double, 2>& angle,
                          double fluidVolume,
                          int timeStepIndex,
                          int reservoirIdx)
{
    // compute the fluid body at the given rotation angles and current fluid volume
    const auto fluidBody = reservoir.computeFluidBody(angle, fluidVolume);

    std::cout << "========================================================================" << std::endl;
    std::cout << "Reservoir " << reservoirIdx << ": optimized reservoirVolume at " << fluidBody.volume << " μl (target: " << fluidVolume << " µl)" << std::endl;
    std::string outputName = "intersections-reservoir_" + std::to_string(reservoirIdx) + '-' + std::to_string(timeStepIndex);

    // for visualization write out the fluid body geometry
    //Dumux::writeIntersections<Point>(fluidBody.triangulation, outputName);

    // the reference points where the inlet/outlet quantities are to be measured
    using Point = Dune::FieldVector<double, 3>;
    static const auto ref0 = Dumux::getParam<Point>("Problem.MeasurementPoint1");
    static const auto ref1 = Dumux::getParam<Point>("Problem.MeasurementPoint2");

    // split volume at y-axis into two reservoirs (approximation of fluid volume that is left for each channel)
    const auto [vol0, vol1] = computeReservoirVolumeWithSplit<Point>(fluidBody.triangulation, 0.5*(ref0[1]+ref1[1]));
    std::cout << "Volume available for Channel (1): " << vol0 << std::endl;
    std::cout << "Volume available for Channel (2): " << vol1 << std::endl;

    // distance of measurement point to water table
    const auto dist0 = std::max(0.0, (fluidBody.waterTablePoint - ref0)*fluidBody.waterTableNormal);
    const auto dist1 = std::max(0.0, (fluidBody.waterTablePoint - ref1)*fluidBody.waterTableNormal);

    // we consider a channel to be "dry" if the water table is below a given threshold
    static const auto dryThreshold = Dumux::getParam<double>("Problem.DryingThresholdInMilliMeter", 1.5);
    const bool dry0 = dist0 < dryThreshold;
    const bool dry1 = dist1 < dryThreshold;

    std::cout << "Water heights:      " << dist0 << ", " << dist1 << " mm "
              << "--> dry?: " << std::boolalpha << dry0 << ", " << dry1 << std::endl;

    // since the platform rotates, the reference points are rotated in space leading to
    // a certain elevation above zero altitude
    const auto rotateX = Dumux::make3DRotation({1.0, 0.0, 0.0}, angle[0]);
    const auto rotateY = Dumux::make3DRotation({0.0, 1.0, 0.0}, angle[1]);
    const auto elevation0 = rotateX(rotateY(ref0))[2];
    const auto elevation1 = rotateX(rotateY(ref1))[2];

    std::cout << "Channel elevations: " << elevation0 << ", " << elevation1 << " mm" << std::endl;

    // The inlet/outlet pressure is computed from the total pressure head (dist + elevation)
    const auto p0 = 1000*9.81*(dist0 + elevation0)*1e-3;
    const auto p1 = 1000*9.81*(dist1 + elevation1)*1e-3;
    std::cout << "Total pressure:     " << p0 << ", " << p1 << " Pa" << std::endl;
    std::cout << "========================================================================" << std::endl;

    {
        std::ofstream metaData(outputName + ".txt");
        metaData << Dumux::Fmt::format("{} {} {} {} {} {} {}\n", angle[0], angle[1], fluidVolume, p0, p1, dry0, dry1);
    }

    struct ReservoirChannels { ChannelState ch0; ChannelState ch1; };
    return ReservoirChannels{ ChannelState{ p0, vol0, dry0 }, ChannelState{ p1, vol1, dry1 } };
}

} // end namespace Dumux::Microfluidic

#endif
