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
#include <tuple>
#include <limits>

#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <dumux/common/timeloop.hh>

#include "reservoir.hh"
#include "channel.hh"
#include "motionfunction.hh"

int main(int argc, char** argv)
{
    // read parameter file "params.input" and command line arguments
    Dumux::Parameters::init(argc, argv);

    // construct a grid representation of the reservoir geometry ("Grid.File")
    Dumux::Microfluidic::Reservoir reservoir{};

    // initial conditions
    const auto reservoirVolume = reservoir.volume();
    const auto initialVolume = Dumux::getParam<double>("Problem.InitialVolumeInMicroLiter", 300.0);
    const auto channelVolume = 0.0; // Consider the channels always filled
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
    const auto channelLength = Dumux::getParam<double>("Problem.ChannelLengthInMeter");
    const auto density = Dumux::getParam<double>("Problem.Density", 1000.0);
    const auto channelCrossSectionArea = Dumux::getParam<double>("Problem.ChannelAreaInSquareMeter");
    const auto zeta = channelTransmissibility*channelLength*density/channelCrossSectionArea;
    std::cout << "Inertial time-scale zeta (s): " << zeta << std::endl;
    std::array<double, 2> fluxInChannelOld{{0.0, 0.0}};
    std::array<double, 2> fluxInChannel{{0.0, 0.0}};

    // max(WSSx)/Q, conversion factor from flow rate to velocity (this is a constant for laminar flow)
    const auto channelWSSFactor = Dumux::getParam<double>("Problem.ChannelWSSFactor");

    // write an output file as text file so that we can make plots over time
    const auto outputFileName = Dumux::getParam<std::string>("Problem.OutputFileName", "output.txt");
    std::ofstream output(outputFileName);
    output << "Time[s] volTotal[μl] volA[μl] volB[μl] flux_ch0[μl/s] flux_ch1[μl/s] max_wss_ch0[Pa] max_wss_ch1[Pa] beta[rad] gamma[rad]\n";

    // parameter for the ad-hoc "capillary stop valve" model
    const double capillaryStopPressure = Dumux::getParam<double>("Problem.CapillaryStopPressure");

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

    // One time step kernel
    auto evolve = [&](const double dt, const double curTime, std::size_t timeStepIndex, const double zeta, bool writeOutput = true)
    {
        // compute current rotation angles of the platform
        double beta = 0.0, gamma = 0.0;
        std::tie(beta, gamma) = motionFunction->rotationAngles(curTime);

        // compute current channel states (pressure, available volume, dry/wet)
        // for the first reservoir use the given angles
        // for the second reservoir we use minus the given angles (geometry mirrored at center of rotation (origin))
        const auto reservoir0State = Dumux::Microfluidic::computeChannelStates(reservoir, { gamma, beta }, volumes[0], timeStepIndex, 0, writeOutput);
        const auto reservoir1State = [&]{
            auto s = Dumux::Microfluidic::computeChannelStates(reservoir, { -gamma, -beta }, volumes[1], timeStepIndex, 1, writeOutput);
            std::swap(s.ch0, s.ch1); // mirrored channel states
            return s;
        }();

        // pressure differences in the two channels based on local elevation and water table height difference
        std::array<double, 2> diffP{{
            reservoir0State.ch0.pressure() - reservoir1State.ch0.pressure(),
            reservoir0State.ch1.pressure() - reservoir1State.ch1.pressure()
        }};

        // Update volumes using pressure gradients and channel resistance (and convert to μl/s)
        // under consideration of inertia
        const auto GLT0 = -1e9*diffP[0]*channelTransmissibility;
        const auto GLT1 = -1e9*diffP[1]*channelTransmissibility;
        fluxInChannel[0] = (fluxInChannelOld[0]*zeta + dt*GLT0)/(zeta + dt);
        fluxInChannel[1] = (fluxInChannelOld[1]*zeta + dt*GLT1)/(zeta + dt);

        ////////////////////////////////////////////////////////////////////////////////
        // Flux limiters ///////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        //
        // Several mechanisms limit the flux that will actually flow:
        //
        // (1) The upstream inlet of the channel has to be wet (connected to the water body)
        //     otherwise the flux is zero.
        //
        // (2) The pressure difference has to be higher than the capillary pressure required
        //     for imbibition of the reservoir (water wants to rather stay in the narrow channel
        //     minimizing liquid-gas interface curvature) so that flow occurs.
        //
        // (3) The maximum amount that can flow in this time step is limited by the available
        //     and connected upstream water body volume
        //
        // (4) The maximum amount that can flow in this time step is limited by the available
        //     space in the downstream reservoir (no overflow)
        //
        //////////////////////////////////////////////////////////////////////////////////

        // flow is only happening if the upstream connection is not dry (1)
        if ((diffP[0] > 0.0 && reservoir0State.ch0.isDry()) || (diffP[0] < 0.0 && reservoir1State.ch0.isDry()))
        {
            std::cout << "Channel (0) stopped by dry upstream inlet" << std::endl;
            fluxInChannel[0] = 0.0;
        }
        // if the downstream is dry flow only happens if the capillary pressure is overcome (piston-imbibition) (2)
        else if ((reservoir1State.ch0.isDry() && diffP[0] > 0.0 && diffP[0] < capillaryStopPressure) || (reservoir0State.ch0.isDry() && diffP[0] < 0.0 && diffP[0] > -capillaryStopPressure))
        {
            std::cout << "Channel (0) stopped by valve pressure: " << capillaryStopPressure << " with ΔP: " << diffP[0] << std::endl;
            fluxInChannel[0] = 0.0;
        }

        // flow is only happening if the upstream connection is not dry (1)
        if ((diffP[1] > 0.0 && reservoir0State.ch1.isDry()) || (diffP[1] < 0.0 && reservoir1State.ch1.isDry()))
        {
            std::cout << "Channel (1) stopped by dry upstream inlet" << std::endl;
            fluxInChannel[1] = 0.0;
        }
        // if the downstream is dry flow only happens if the capillary pressure is overcome (piston-imbibition) (2)
        else if ((reservoir1State.ch1.isDry() && diffP[1] > 0.0 && diffP[1] < capillaryStopPressure) || (reservoir0State.ch1.isDry() && diffP[1] < 0.0 && diffP[1] > -capillaryStopPressure))
        {
            std::cout << "Channel (1) stopped by valve pressure: " << capillaryStopPressure << " with ΔP: " << diffP[1] << std::endl;
            fluxInChannel[1] = 0.0;
        }

        // flow from reservoir 0 to reservoir 1 in channel 0 (3)
        constexpr double fluxEps = 1e-5;
        if (fluxInChannel[0] < -fluxEps)
            std::cout << "Volume available upstream for Channel (0): " << reservoir0State.ch0.availableFluidVolume() << std::endl;
        else if (fluxInChannel[0] > fluxEps)
            std::cout << "Volume available upstream for Channel (0): " << reservoir1State.ch0.availableFluidVolume() << std::endl;

        if (fluxInChannel[1] < -fluxEps)
            std::cout << "Volume available upstream for Channel (1): " << reservoir0State.ch1.availableFluidVolume() << std::endl;
        else if (fluxInChannel[1] > fluxEps)
            std::cout << "Volume available upstream for Channel (1): " << reservoir1State.ch1.availableFluidVolume() << std::endl;

        if (fluxInChannel[0] < -fluxEps && fluxInChannel[0] < -reservoir0State.ch0.availableFluidVolume()/dt)
        {
            fluxInChannel[0] = -reservoir0State.ch0.availableFluidVolume()/dt;
            std::cout << "Channel (0) flux limited because too little upstream water is available" << std::endl;
        }
        // flow from reservoir 1 to reservoir 0 in channel 0 (3)
        else if (fluxInChannel[0] > fluxEps && fluxInChannel[0] > reservoir1State.ch0.availableFluidVolume()/dt)
        {
            fluxInChannel[0] = reservoir1State.ch0.availableFluidVolume()/dt;
            std::cout << "Channel (0) flux limited because too little upstream water is available" << std::endl;
        }

        // flow from reservoir 0 to reservoir 1 in channel 1 (3)
        if (fluxInChannel[1] < -fluxEps && fluxInChannel[1] < -reservoir0State.ch1.availableFluidVolume()/dt)
        {
            fluxInChannel[1] = -reservoir0State.ch1.availableFluidVolume()/dt;
            std::cout << "Channel (1) flux limited because too little upstream water is available" << std::endl;
        }
        // flow from reservoir 1 to reservoir 0 in channel 1 (3)
        else if (fluxInChannel[1] > fluxEps && fluxInChannel[1] > reservoir1State.ch1.availableFluidVolume()/dt)
        {
            fluxInChannel[1] = reservoir1State.ch1.availableFluidVolume()/dt;
            std::cout << "Channel (1) flux limited because too little upstream water is available" << std::endl;
        }

        const auto netFluxPredict = fluxInChannel[0] + fluxInChannel[1];
        // make sure only as much flows as is actually there; (3) and (4)
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

        if (fluxInChannel[0] < -1e-7)
            std::cout << "\x1b[31m" << "Channel (0) flow R0 " << std::string(int(std::abs(fluxInChannel[0])), '>') << "> R1 (ΔP: " << diffP[0] << " Pa, Q: " << fluxInChannel[0] <<  " µl/s, V: " << std::abs(fluxInChannel[0]*dt) <<  " µl)" << "\033[0m" << "\n";
        else if (fluxInChannel[0] > 1e-7)
            std::cout << "\x1b[31m" << "Channel (0) flow R0 <" << std::string(int(std::abs(fluxInChannel[0])), '<') << " R1 (ΔP: " << diffP[0] << " Pa, Q: " << fluxInChannel[0] <<  " µl/s, V: " << std::abs(fluxInChannel[0]*dt) <<  " µl)" << "\033[0m" << "\n";
        else
            std::cout << "\x1b[31m" << "Channel (0) no flow" << "\033[0m" << "\n";

        if (fluxInChannel[1] < -1e-7)
            std::cout << "\x1b[31m" << "Channel (1) flow R0 " << std::string(int(std::abs(fluxInChannel[1])), '>') << "> R1 (ΔP: " << diffP[1] << " Pa, Q: " << fluxInChannel[1] <<  " µl/s, V: " << std::abs(fluxInChannel[1]*dt) <<  " µl)" << "\033[0m" << "\n";
        else if (fluxInChannel[1] > 1e-7)
            std::cout << "\x1b[31m" << "Channel (1) flow R0 <" << std::string(int(std::abs(fluxInChannel[1])), '<') << " R1 (ΔP: " << diffP[1] << " Pa, Q: " << fluxInChannel[1] <<  " µl/s, V: " << std::abs(fluxInChannel[1]*dt) <<  " µl)" << "\033[0m" << "\n";
        else
            std::cout << "\x1b[31m" << "Channel (1) no flow" << "\033[0m" << "\n";

        volumes[0] += netVolumeExchange;
        volumes[1] -= netVolumeExchange;

        std::cout << "Time: " << curTime << std::endl;

        if (writeOutput)
        {
            output << curTime << " "
                << volumes[0] + volumes[1] + 2*channelVolume << " "
                << volumes[0] << " "
                << volumes[1] << " "
                << fluxInChannel[0] << " "
                << fluxInChannel[1] << " "
                << std::abs(fluxInChannel[0]*1e-9*channelWSSFactor) << " "
                << std::abs(fluxInChannel[1]*1e-9*channelWSSFactor) << " "
                << beta << " "
                << gamma
                << std::endl;
        }

        // update for inertia model so we can compute the rate of change of the flux
        fluxInChannelOld = fluxInChannel;
    };


    // Initialization phase > Simulate until water is in equilibrium and the flux is zero
    auto oldVolumes = volumes;
    while (true)
    {
        const auto dt = timeLoop->timeStepSize();

        // stay at starting time and neglect inertial effects (we want to get there fast)
        evolve(dt, 0.0, 0, 0.0, false);

        if (std::abs(oldVolumes[0] - volumes[0]) < 0.01)
            break;

        oldVolumes = volumes;
    }

    // Actual simulation
    timeLoop->start(); do
    {
        const auto curTime = timeLoop->time();
        const auto timeStepIndex = timeLoop->timeStepIndex();
        const auto dt = timeLoop->timeStepSize();

        // do one time step
        evolve(dt, curTime, timeStepIndex, zeta, true);

        // go to next time step
        timeLoop->advanceTimeStep();

        // report statistics of this time step
        timeLoop->reportTimeStep();
        std::cout << std::endl;

    } while (!timeLoop->finished());

    timeLoop->finalize();

    return 0;
}
