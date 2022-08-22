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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/

#include <config.h>

#include <vector>
#include <array>

#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <dumux/common/timeloop.hh>

#include "reservoir.hh"
#include "channel.hh"
#include "motionfunction.hh"

std::array<double, 2> pressureGradients(const Dumux::Microfluidic::Reservoir& reservoir,
                                        const std::array<double, 2> angle,
                                        const std::array<double, 2> volumes,
                                        int timeStepIndex)
{
    // for the first reservoir use the given angles
    const auto reservoir0State = Dumux::Microfluidic::computeChannelStates(reservoir, angle, volumes[0], timeStepIndex, 0);
    // for the second reservoir we use minus the given angles (geometry mirrored at center of rotation (origin))
    const auto reservoir1State = [&]{
        auto s = Dumux::Microfluidic::computeChannelStates(reservoir, { -angle[0], -angle[1] }, volumes[1], timeStepIndex, 1);
        std::swap(s.ch0, s.ch1); // mirrored channel states
        return s;
    }();

    auto diffP0 = reservoir0State.ch0.pressure() - reservoir1State.ch0.pressure();
    auto diffP1 = reservoir0State.ch1.pressure() - reservoir1State.ch1.pressure();

    // parameter for the ad-hoc "capillary stop valve" model
    static const double capillaryStopPressure = Dumux::getParam<double>("Problem.CapillaryStopPressure");

    // flow is only happening if the upstream connection is not dry
    if ((diffP0 > 0.0 && reservoir0State.ch0.isDry()) || (diffP0 < 0.0 && reservoir1State.ch0.isDry()))
    {
        std::cout << "Channel (1) stopped by dry upstream inlet" << std::endl;
        diffP0 = 0.0;
    }
    // if the downstream is dry flow only happens if the capillary pressure is overcome (piston-imbibition)
    else if ((reservoir1State.ch0.isDry() && diffP0 > 0.0 && diffP0 < capillaryStopPressure) || (reservoir0State.ch0.isDry() && diffP0 < 0.0 && diffP0 > -capillaryStopPressure))
    {
        std::cout << "Channel (1) stopped by valve pressure: " << capillaryStopPressure << " with ΔP: " << diffP0 << std::endl;
        diffP0 = 0.0;
    }
    // flow is only happening if the upstream connection is not dry
    if ((diffP1 > 0.0 && reservoir0State.ch1.isDry()) || (diffP1 < 0.0 && reservoir1State.ch1.isDry()))
    {
        std::cout << "Channel (2) stopped by dry upstream inlet" << std::endl;
        diffP1 = 0.0;
    }
    // if the downstream is dry flow only happens if the capillary pressure is overcome (piston-imbibition)
    else if ((reservoir1State.ch1.isDry() && diffP1 > 0.0 && diffP1 < capillaryStopPressure) || (reservoir0State.ch1.isDry() && diffP1 < 0.0 && diffP1 > -capillaryStopPressure))
    {
        std::cout << "Channel (2) stopped by valve pressure: " << capillaryStopPressure << " with ΔP: " << diffP1 << std::endl;
        diffP1 = 0.0;
    }

    return { diffP0, diffP1 };
}

int main(int argc, char** argv)
{
    // read parameter file "params.input" and command line arguments
    Dumux::Parameters::init(argc, argv);

    // construct a grid representation of the reservoir geometry
    Dumux::Microfluidic::Reservoir reservoir{};

    // initial conditions
    const auto reservoirVolume = reservoir.volume();
    const auto initialVolume = Dumux::getParam<double>("Problem.InitialVolumeInMicroLiter", 300.0);
    const auto channelVolume = Dumux::getParam<double>("Problem.SingleChannelVolumeInMicroLiter", 16.735);
    std::cout << "Initial reservoirVolume (µl): " << initialVolume << ", channel reservoirVolume (µl): " << channelVolume << std::endl;
    if (auto maxV = 2*channelVolume + 2*reservoirVolume; initialVolume > maxV)
        DUNE_THROW(Dune::IOError, "The maximum volume this chip can fit is " << maxV << " µl");

    // we start with as much as we can fit into the first volume and put the rest in the second volume
    const auto firstReservoirVolume = std::min(reservoirVolume, initialVolume - 2*channelVolume);
    const auto secondReservoirVolume = std::max(0.0, initialVolume - 2*channelVolume - firstReservoirVolume);
    // here we store the current primary variables (volume in each reservoir, V1 and V2)
    std::array<double, 2> volumes({ firstReservoirVolume, secondReservoirVolume });
    std::cout << "Available reservoirVolume (µl): " << reservoirVolume << std::endl;
    std::cout << "Initial volumes (µl): R0: " << volumes[0] << ", R1: " << volumes[1] << " µl" << std::endl;

    // m^3/(s*Pa)
    const auto channelTransmissibility = Dumux::getParam<double>("Problem.ChannelTransmissibility");
    // m/(s*Pa), conversion factor from flow rate to velocity (this is a constant for Stokes flow)
    const auto channelVelFactor = Dumux::getParam<double>("Problem.ChannelVelFactor");

    // write an output file as text file so that we can make plots over time
    const auto outputFileName = Dumux::getParam<std::string>("Problem.OutputFileName", "output.txt");
    std::ofstream output(outputFileName);
    output << "Time[s] volTotal[μl] volA[μl] volB[μl] flux_ch0[μl/s] flux_ch1[μl/s] maxv_ch0[m/s] maxv_ch1[m/s] beta[rad] gamma[rad]\n";

    // parameters for ad-hoc model of inertia effects
    // maximum flux change per second
    const double maxFluxChangePerSecond = Dumux::getParam<double>("Problem.MaxFluxChangePerSecond", 10);
    std::array<double, 2> fluxInChannelOld{{0.0, 0.0}};
    std::array<double, 2> fluxInChannel{{0.0, 0.0}};
    std::array<double, 2> fluxInChannelDeriv{{0.0, 0.0}};

    // the motion function (using constant tilt here / can be replaced by something else in the future)
    const auto rotationsPerSecond = Dumux::getParam<double>("Problem.RotationsPerMinute", 4.0)/60.0;
    std::shared_ptr<Dumux::Microfluidic::MotionFunction> motionFunction
        = std::make_unique<Dumux::Microfluidic::ConstantTiltMotionFunction>(rotationsPerSecond);

    // Read parameters and setup the time loop
    const auto cycles = Dumux::getParam<double>("TimeLoop.Cycles", 1.0);
    const auto tEnd = Dumux::getParam<double>("TimeLoop.TEnd", cycles/rotationsPerSecond);
    const auto stepsPerCycle = Dumux::getParam<double>("TimeLoop.StepsPerCycle", 100);
    const auto dt = Dumux::getParam<double>("TimeLoop.Dt", tEnd/cycles/stepsPerCycle);
    auto timeLoop = std::make_shared<Dumux::TimeLoop<double>>(0.0, dt, tEnd);

    timeLoop->start(); do
    {
        // compute current rotation angles of the platform
        const auto curTime = timeLoop->time();
        const auto [beta, gamma] = motionFunction->rotationAngles(curTime);

        // compute current pressure gradients for the channels
        const auto diffP = pressureGradients(reservoir, { gamma, beta }, volumes, timeLoop->timeStepIndex());

        // Update volumes using pressure gradients and channel resistance
        // convert to μl/s
        fluxInChannel[0] = -1e9*diffP[0]*channelTransmissibility;
        fluxInChannel[1] = -1e9*diffP[1]*channelTransmissibility;
        const auto dt = timeLoop->timeStepSize();

        // check that the flux doesn't grow too much per time (inertia)
        {
            fluxInChannelDeriv[0] = (fluxInChannel[0] - fluxInChannelOld[0])/dt;
            fluxInChannelDeriv[1] = (fluxInChannel[1] - fluxInChannelOld[1])/dt;
            const auto regularizeFlux = [&](int i){
                if (fluxInChannel[i] > 0 && fluxInChannelDeriv[i] > 0)
                {
                    if (fluxInChannelDeriv[i] > maxFluxChangePerSecond)
                    {
                        fluxInChannel[i] = fluxInChannelOld[i] + dt*maxFluxChangePerSecond;
                        return true;
                    }
                }
                else if (fluxInChannel[i] < 0 && fluxInChannelDeriv[i] < 0)
                {
                    if (fluxInChannelDeriv[i] < -maxFluxChangePerSecond)
                    {
                        fluxInChannel[i] = fluxInChannelOld[i] - dt*maxFluxChangePerSecond;
                        return true;
                    }
                }
                return false;
            };

            bool regularized = regularizeFlux(0);
            if (regularized)
                std::cout << "Regularized Channel (1) flux/dt: " << fluxInChannelDeriv[0] << " max: " << maxFluxChangePerSecond << " µl/s^2" << std::endl;
            regularized = regularizeFlux(1);
            if (regularized)
                std::cout << "Regularized Channel (2) flux/dt: " << fluxInChannelDeriv[1] << " max: " << maxFluxChangePerSecond << " µl/s^2" << std::endl;
        }

        const auto netFluxPredict = fluxInChannel[0] + fluxInChannel[1];
        // make sure only as much flows as is actually there
        const auto predictedNetVolumeExchange = netFluxPredict*timeLoop->timeStepSize();
        const auto netVolumeExchange = std::clamp(predictedNetVolumeExchange,
            // minimum is constrained by what's left in reservoirVolume 0 or can be added to reservoirVolume 1
            std::max(-volumes[0], volumes[1]-reservoirVolume),
            // maximum is constrained by what's left in reservoirVolume 1 or can be added to reservoirVolume 0
            std::min(volumes[1], reservoirVolume-volumes[0])
        );

        if (predictedNetVolumeExchange < std::max(-volumes[0], volumes[1]-reservoirVolume))
        {
            std::cout << "\x1b[31m" << "Restricted net flux because it was too negative (w.r.t. vol0)" << "\033[0m" << "\n";
            std::cout << "\x1b[31m" << predictedNetVolumeExchange << "< max( " << -volumes[0] << ", " << volumes[1]-reservoirVolume << ")\033[0m" << "\n";
            std::cout << "\x1b[31m" << predictedNetVolumeExchange << "< " << std::max(-volumes[0], volumes[1]-reservoirVolume) << "\033[0m" << "\n";
            const auto fraction = std::clamp(netVolumeExchange/predictedNetVolumeExchange, 0.0, 1.0);
            fluxInChannel[0] *= fraction;
            fluxInChannel[1] *= fraction;
        }

        if (predictedNetVolumeExchange > std::min(volumes[1], reservoirVolume-volumes[0]))
        {
            std::cout << "\x1b[31m" << "Restricted net flux because it was too positive (w.r.t. vol0)" << "\033[0m" << "\n";
            std::cout << "\x1b[31m" << predictedNetVolumeExchange << "> min( " << volumes[1] << ", " << reservoirVolume-volumes[0] << ")\033[0m" << "\n";
            std::cout << "\x1b[31m" << predictedNetVolumeExchange << "> " << std::min(volumes[1], reservoirVolume-volumes[0]) << "\033[0m" << "\n";
            const auto fraction = std::clamp(netVolumeExchange/predictedNetVolumeExchange, 0.0, 1.0);
            fluxInChannel[0] *= fraction;
            fluxInChannel[1] *= fraction;
        }

        if (diffP[0] > 1e-7)
            std::cout << "\x1b[31m" << "Channel (1) flow R0 " << std::string(int(std::abs(fluxInChannel[0])), '>') << "> R1 (ΔP: " << diffP[0] << " Pa, Q: " << fluxInChannel[0] <<  " µl/s)" << "\033[0m" << "\n";
        else if (diffP[0] < -1e-7)
            std::cout << "\x1b[31m" << "Channel (1) flow R0 <" << std::string(int(std::abs(fluxInChannel[0])), '<') << " R1 (ΔP: " << diffP[0] << " Pa, Q: " << fluxInChannel[0] <<  " µl/s)" << "\033[0m" << "\n";
        else
            std::cout << "\x1b[31m" << "Channel (1) no flow" << "\033[0m" << "\n";

        if (diffP[1] > 1e-7)
            std::cout << "\x1b[31m" << "Channel (2) flow R0 " << std::string(int(std::abs(fluxInChannel[1])), '>') << "> R1 (ΔP: " << diffP[1] << " Pa, Q: " << fluxInChannel[1] <<  " µl/s)" << "\033[0m" << "\n";
        else if (diffP[1] < -1e-7)
            std::cout << "\x1b[31m" << "Channel (2) flow R0 <" << std::string(int(std::abs(fluxInChannel[1])), '<') << " R1 (ΔP: " << diffP[1] << " Pa, Q: " << fluxInChannel[1] <<  " µl/s)" << "\033[0m" << "\n";
        else
            std::cout << "\x1b[31m" << "Channel (2) no flow" << "\033[0m" << "\n";

        volumes[0] += netVolumeExchange;
        volumes[1] -= netVolumeExchange;

        std::cout << "Time: " << curTime << std::endl;
        output << curTime << " "
               << volumes[0] + volumes[1] + 2*channelVolume << " "
               << volumes[0] << " "
               << volumes[1] << " "
               << fluxInChannel[0] << " "
               << fluxInChannel[1] << " "
               << -(fluxInChannel[0]*1e-9/channelTransmissibility*channelVelFactor) << " "
               << -(fluxInChannel[1]*1e-9/channelTransmissibility*channelVelFactor) << " "
               << beta << " "
               << gamma
               << std::endl;

        // update for inertia model so we can compute the rate of change of the flux
        fluxInChannelOld = fluxInChannel;

        // go to next time step
        timeLoop->advanceTimeStep();

        // report statistics of this time step
        timeLoop->reportTimeStep();
        std::cout << std::endl;

    } while (!timeLoop->finished());

    timeLoop->finalize();

    return 0;
}
