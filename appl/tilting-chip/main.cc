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

#include <dune/common/timer.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/common/math.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/nonlinear/findscalarroot.hh>

#include <dumux/io/grid/gridmanager_alu.hh>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>

#include <test/geometry/writetriangulation.hh>

#include "intersections.hh"
#include "rotation.hh"
#include "reservoir.hh"
#include "volume.hh"

// compute the pressure at channel connections of a reservoir
template<class Point, class Tree>
auto computeChannelStates(const Tree& aabbTree, std::vector<Point> corners, std::array<double, 2> angle, double reservoirVolume, int timeStepIndex, int reservoirIdx)
{
    // adjust the top plane of the helper hexahedron to simulate the water table angle
    Point rotationCenter(0.0);
    for (int i = 4; i < 8; ++i)
        rotationCenter += corners[i];
    rotationCenter /= 4.0;

    for (int i = 4; i < 8; ++i)
    {
        const auto heightDiffX = std::tan(angle[0])*(corners[i][1] - rotationCenter[1]);
        const auto heightDiffY = std::tan(angle[1])*(corners[i][0] - rotationCenter[0]);
        corners[i][2] += heightDiffX + heightDiffY;
    }

    using HexGeometry = Dune::MultiLinearGeometry<double, 3, 3>;
    auto hexCut = HexGeometry(Dune::GeometryTypes::cube(3), corners);

    // compute intersections of the helper hexahedron with the reservoir geometry
    auto treeIntersections = Dumux::intersectingEntities(hexCut, aabbTree);
    auto intersections = Dumux::convertIntersections<Point>(treeIntersections);

    // find the water table height by optimzing the reservoirVolume to the given one
    auto localCorners = corners; // avoid reallocation every iteration
    const auto residual = [&](const double h)
    {
        // shift water table height
        localCorners = corners;
        for (int i = 4; i < 8; ++i)
            localCorners[i][2] += h;

        hexCut = HexGeometry(Dune::GeometryTypes::cube(3), localCorners);
        Dumux::intersectingEntities(treeIntersections, hexCut, aabbTree);
        Dumux::convertIntersections<Point>(intersections, treeIntersections);
        const auto iVol = Dumux::computeReservoirVolume<Point>(intersections);
        // std::cout << "h: " << h << "-> vol: " << iVol << " / residual: " << iVol - reservoirVolume << std::endl;
        return iVol - reservoirVolume;
    };
    double h = 0.0;
    try {
        h = Dumux::findScalarRootBrent(-11.0, 5.0, residual, 1e-2, 2000);
    }
    catch (const Dune::InvalidStateException& e)
    {
        h = Dumux::findScalarRootBrent(-11.1, 5.1, residual, 1e-2, 2000);
    }
    auto vol = residual(h) + reservoirVolume;

    std::cout << "========================================================================" << std::endl;
    std::cout << "Reservoir " << reservoirIdx << ": optimized reservoirVolume at " << vol << " μl (target: " << reservoirVolume << " µl)" << std::endl;
    std::string outputName = "intersections-reservoir_" + std::to_string(reservoirIdx) + '-' + std::to_string(timeStepIndex);
    //writeIntersections<Point>(intersections, outputName);

    static const auto ref0 = Dumux::getParam<Point>("Problem.MeasurementPoint1");
    static const auto ref1 = Dumux::getParam<Point>("Problem.MeasurementPoint2");

    // split volume at y-axis into two reservoirs (approximation of reservoir volume that is left)
    const auto [vol0, vol1] = Dumux::computeReservoirVolumeWithSplit<Point>(intersections, 0.5*(ref0[1]+ref1[1]));
    std::cout << "Volume available for channel 0: " << vol0 << std::endl;
    std::cout << "Volume available for channel 1: " << vol1 << std::endl;

    localCorners = corners;
    for (int i = 4; i < 8; ++i)
        localCorners[i][2] += h;

    const auto ab = localCorners[5] - localCorners[4];
    const auto ac = localCorners[6] - localCorners[4];
    auto normal = Dumux::crossProduct(ab, ac);
    normal /= normal.two_norm();

    const auto dist0 = std::max(0.0, (localCorners[4] - ref0)*normal);
    const auto dist1 = std::max(0.0, (localCorners[4] - ref1)*normal);

    static const auto dryThreshold = Dumux::getParam<double>("Problem.DryingThresholdInMilliMeter", 1.5);
    const bool dry0 = dist0 < dryThreshold;
    const bool dry1 = dist1 < dryThreshold;

    std::cout << "Water heights:      " << dist0 << ", " << dist1 << " mm "
              << "--> dry?: " << std::boolalpha << dry0 << ", " << dry1 << std::endl;

    const auto rotateX = Dumux::make3DRotation({1.0, 0.0, 0.0}, angle[0]);
    const auto rotateY = Dumux::make3DRotation({0.0, 1.0, 0.0}, angle[1]);
    const auto elevation0 = rotateX(rotateY(ref0))[2];
    const auto elevation1 = rotateX(rotateY(ref1))[2];

    std::cout << "Channel elevations: " << elevation0 << ", " << elevation1 << " mm" << std::endl;

    const auto p0 = 1000*9.81*(dist0 + elevation0)*1e-3;
    const auto p1 = 1000*9.81*(dist1 + elevation1)*1e-3;
    std::cout << "Total pressure:     " << p0 << ", " << p1 << " Pa" << std::endl;
    std::cout << "========================================================================" << std::endl;

    {
        std::ofstream metaData(outputName + ".txt");
        metaData << Dumux::Fmt::format("{} {} {} {} {} {} {}\n", angle[0], angle[1], reservoirVolume, p0, p1, dry0, dry1);
    }

    struct ChannelState { double pressure; bool dry; double fractionalVolume; };
    struct ReservoirChannels { ChannelState ch0; ChannelState ch1; };

    return ReservoirChannels{ ChannelState{ p0, dry0, vol0 }, ChannelState{ p1, dry1, vol1 } };
}

template<class Point, class Tree>
std::array<double, 2> pressureGradients(const Tree& aabbTree, const std::vector<Point>& corners, std::array<double, 2> angle, std::array<double, 2> volumes, int timeStepIndex)
{
    static const double capillaryStopPressure = Dumux::getParam<double>("Problem.CapillaryStopPressure");
    // for the first reservoir use the given angles
    const auto reservoir0State = computeChannelStates(aabbTree, corners, angle, volumes[0], timeStepIndex, 0);
    // for the second reservoir we use minus the given angles (geometry mirrored at center of rotation (origin))
    const auto reservoir1State = [&]{
        auto s = computeChannelStates(aabbTree, corners, { -angle[0], -angle[1] }, volumes[1], timeStepIndex, 1);
        std::swap(s.ch0, s.ch1); // mirrored channel states
        return s;
    }();

    auto diffP0 = reservoir0State.ch0.pressure - reservoir1State.ch0.pressure;
    auto diffP1 = reservoir0State.ch1.pressure - reservoir1State.ch1.pressure;

    // flow is only happening if the upstream connection is not dry
    if ((diffP0 > 0.0 && reservoir0State.ch0.dry) || (diffP0 < 0.0 && reservoir1State.ch0.dry))
    {
        std::cout << "Channel (1) stopped by dry upstream inlet" << std::endl;
        diffP0 = 0.0;
    }
    // if the downstream is dry flow only happens if the capillary pressure is overcome (piston-imbibition)
    else if ((reservoir1State.ch0.dry && diffP0 > 0.0 && diffP0 < capillaryStopPressure) || (reservoir0State.ch0.dry && diffP0 < 0.0 && diffP0 > -capillaryStopPressure))
    {
        std::cout << "Channel (1) stopped by valve pressure: " << capillaryStopPressure << " with ΔP: " << diffP0 << std::endl;
        diffP0 = 0.0;
    }
    // flow is only happening if the upstream connection is not dry
    if ((diffP1 > 0.0 && reservoir0State.ch1.dry) || (diffP1 < 0.0 && reservoir1State.ch1.dry))
    {
        std::cout << "Channel (2) stopped by dry upstream inlet" << std::endl;
        diffP1 = 0.0;
    }
    // if the downstream is dry flow only happens if the capillary pressure is overcome (piston-imbibition)
    else if ((reservoir1State.ch1.dry && diffP1 > 0.0 && diffP1 < capillaryStopPressure) || (reservoir0State.ch1.dry && diffP1 < 0.0 && diffP1 > -capillaryStopPressure))
    {
        std::cout << "Channel (2) stopped by valve pressure: " << capillaryStopPressure << " with ΔP: " << diffP1 << std::endl;
        diffP1 = 0.0;
    }

    return { diffP0, diffP1 };
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // read parameter file "params.input" and command line arguments
    Parameters::init(argc, argv);

    // construct a grid representation of the reservoir geometry
    Microfluidic::Reservoir reservoir{};
    auto [lowerLeft, upperRight] = reservoir.boundingBox();

    const auto height = upperRight[2]-lowerLeft[2];
    lowerLeft[2] -= height;
    upperRight[2] += height;

    using Point = Dune::FieldVector<double, 3>;
    std::vector<Point> corners{
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
        corners[i][2] -= 1.0*height;

    // rotate once around the clock at a given angle and compute the water body geometry for constant reservoirVolume
    const auto angleDegree = getParam<double>("Problem.Angle", 18.0);
    const auto theta = angleDegree/180.0*M_PI;
    const auto sinTheta = std::sin(theta);

    const auto rotationsPerSecond = getParam<double>("Problem.RotationsPerMinute", 4.0)/60.0;
    const auto cycles = getParam<double>("TimeLoop.Cycles", 1.0);
    const auto tEnd = getParam<double>("TimeLoop.TEnd", cycles/rotationsPerSecond);
    const auto stepsPerCycle = getParam<double>("TimeLoop.StepsPerCycle", 100);
    const auto dt = getParam<double>("TimeLoop.Dt", tEnd/cycles/stepsPerCycle);
    auto timeLoop = std::make_shared<TimeLoop<double>>(0.0, dt, tEnd);

    const auto reservoirVolume = reservoir.volume();
    const auto initialVolume = getParam<double>("Problem.InitialVolumeInMicroLiter", 300.0);
    const auto channelVolume = getParam<double>("Problem.SingleChannelVolumeInMicroLiter", 16.735);
    std::cout << "Initial reservoirVolume (µl): " << initialVolume << ", channel reservoirVolume (µl): " << channelVolume << std::endl;
    if (auto maxV = 2*channelVolume + 2*reservoirVolume; initialVolume > maxV)
        DUNE_THROW(Dune::IOError, "The maximum volume this chip can fit is " << maxV << " µl");
    // we start with as much as we can fit into the first volume and put the rest in the second volume
    const auto firstReservoirVolume = std::min(reservoirVolume, initialVolume - 2*channelVolume);
    const auto secondReservoirVolume = std::max(0.0, initialVolume - 2*channelVolume - firstReservoirVolume);
    std::array<double, 2> volumes({ firstReservoirVolume, secondReservoirVolume });
    std::cout << "Available reservoirVolume (µl): " << reservoirVolume << std::endl;
    std::cout << "Initial volumes (µl): R0: " << volumes[0] << ", R1: " << volumes[1] << " µl" << std::endl;

    // m^3/(s*Pa)
    const auto channelTransmissibility = getParam<double>("Problem.ChannelTransmissibility", 9e-10);
    // m/(s*Pa)
    const auto channelVelFactor = getParam<double>("Problem.ChannelVelFactor", 0.026e-1);
    // maximum flux change per second
    const double maxFluxChangePerSecond = getParam<double>("Problem.MaxFluxChangePerSecond", 10);

    const auto outputFileName = getParam<std::string>("Problem.OutputFileName", "output.txt");
    std::ofstream output(outputFileName);
    output << "Time[s] volTotal[μl] volA[μl] volB[μl] flux_ch0[μl/s] flux_ch1[μl/s] maxv_ch0[m/s] maxv_ch1[m/s] beta[rad] gamma[rad]\n";
    std::array<double, 2> fluxInChannelOld{{0.0, 0.0}};
    std::array<double, 2> fluxInChannel{{0.0, 0.0}};
    std::array<double, 2> fluxInChannelDeriv{{0.0, 0.0}};
    timeLoop->start(); do
    {
        // compute current angles
        const auto curTime = timeLoop->time();
        const auto t = curTime*rotationsPerSecond*2.0*M_PI;
        const auto sinT = std::sin(t);
        const auto cosT = std::cos(t);
        const auto gamma = std::asin(-sinTheta*sinT);
        const auto beta = std::asin(-sinTheta*cosT/std::sqrt(1.0 - sinTheta*sinTheta*sinT*sinT));

        // compute current pressure gradients
        // regularized reservoirVolume so that the bracket algorithm for finding water heights is guaranteed to work
        auto volRegu = volumes;
        volRegu[0] = std::clamp(volRegu[0], 1e-2, reservoirVolume-1e-2);
        volRegu[1] = std::clamp(volRegu[1], 1e-2, reservoirVolume-1e-2);
        const auto diffP = pressureGradients(reservoir.boundingBoxTree(), corners, { gamma, beta }, volRegu, timeLoop->timeStepIndex());

        // Update volumes using pressure gradients and channel resistance
        // convert to μl/s
        fluxInChannel[0] = -1e9*diffP[0]*channelTransmissibility;
        fluxInChannel[1] = -1e9*diffP[1]*channelTransmissibility;
        const auto dt = timeLoop->timeStepSize();

        // check that the flux doesn't grow too much per time (inertia)
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

        std::cout << "Time: " << curTime << ", t: " << t << std::endl;
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