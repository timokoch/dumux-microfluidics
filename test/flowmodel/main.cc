
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
#include <dumux/geometry/intersectingentities.hh>

#include <test/geometry/writetriangulation.hh>

template<class Point, class TreeIntersections>
auto convertIntersections(const TreeIntersections& treeIntersections)
{
    std::vector<std::vector<Point>> intersections;
    intersections.reserve(treeIntersections.size());
    for (const auto& is : treeIntersections)
        intersections.emplace_back(std::vector<Point>(is.corners()));
    return intersections;
}

template<class Point>
void writeIntersections(const std::vector<std::vector<Point>>& intersections, const std::string& name)
{
    // Dune::Timer timer;
    // std::cout << "Writing " << intersections.size() << " intersections to file ...";
    Dumux::writeVTUTetrahedron(intersections, name);
    // std::cout << " done ( " << timer.elapsed() << " seconds)." << std::endl;
}

template<class Point>
void writeIntersections(const std::vector<std::array<Point, 4>>& intersections, const std::string& name)
{
    // Dune::Timer timer;
    // std::cout << "Writing " << intersections.size() << " intersections to file ...";
    Dumux::writeVTUTetrahedron(intersections, name);
    // std::cout << " done ( " << timer.elapsed() << " seconds)." << std::endl;
}

template<class Point>
double computeVolume(const std::vector<std::vector<Point>>& intersections)
{
    using Simplex3Geometry = Dune::AffineGeometry<double, 3, 3>;
    double volume = 0.0;
    for (const auto& tet : intersections)
    {
        const auto tetGeo = Simplex3Geometry(
            Dune::GeometryTypes::simplex(3), std::array<Point, 4>{{ tet[0], tet[1], tet[2], tet[3] }}
        );
        volume += tetGeo.volume();
    }
    return volume;
}

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

template<class Geometry>
void writeHex(const Geometry& geometry, const std::string& filename)
{
    std::ofstream fout(filename + ".vtu");
    fout << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n"
         << "  <UnstructuredGrid>\n"
         << "    <Piece NumberOfPoints=\"" << 8 << "\" NumberOfCells=\"" << 1 << "\">\n"
         << "      <Points>\n"
         << "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    constexpr std::array<int, 8> map{{ 0, 1, 3, 2, 4, 5, 7, 6 }};
    for (int i = 0; i < geometry.corners(); ++i)
        fout << geometry.corner(map[i]) << " ";

    fout << '\n';
    fout << "        </DataArray>\n"
         << "      </Points>\n"
         << "      <Cells>\n"
         << "        <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    fout << "        0 1 2 3 4 5 6 7\n";
    fout << "        </DataArray>\n";
    fout << "        <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    fout << "        8\n";
    fout << "        </DataArray>\n";
    fout << "        <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    fout << "        12\n";
    fout << "        </DataArray>\n"
         << "      </Cells>\n"
         << "    </Piece>\n"
         << "</UnstructuredGrid>\n"
         << "</VTKFile>\n";
}

// compute the pressure at channel connections of a reservoir
template<class Point, class Tree>
auto computeChannelStates(const Tree& aabbTree, std::vector<Point> corners, std::array<double, 2> angle, double volume, int timeStepIndex, int reservoirIdx)
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
    auto treeIntersections = intersectingEntities(hexCut, aabbTree);
    auto intersections = convertIntersections<Point>(treeIntersections);

    // find the water table height by optimzing the volume to the given one
    const auto residual = [&](const double h)
    {
        // shift water table height
        auto localCorners = corners;
        for (int i = 4; i < 8; ++i)
            localCorners[i][2] += h;

        hexCut = HexGeometry(Dune::GeometryTypes::cube(3), localCorners);
        treeIntersections = intersectingEntities(hexCut, aabbTree);
        intersections = convertIntersections<Point>(treeIntersections);
        const auto iVol = computeVolume<Point>(intersections);
        //std::cout << "h: " << h << "-> vol: " << iVol << " / residual: " << iVol - volume << std::endl;
        return iVol - volume;
    };

    auto h = Dumux::findScalarRootBrent(-11.0, 0.0, residual, 1e-2);
    auto vol = residual(h) + volume;

    std::cout << "========================================================================" << std::endl;
    std::cout << "Reservoir " << reservoirIdx << ": optimized volume at " << vol << " μl" << std::endl;
    std::string outputName = "intersections-reservoir_" + std::to_string(reservoirIdx) + '-' + std::to_string(timeStepIndex);
    writeIntersections<Point>(intersections, outputName);

    Point ref0({ 0.5*(3.462435334542 + 3.798170069313) + 7.3, 0.5*(-6.14367530233 -7.067145485502) - 2.0, 0.0 });
    Point ref1({ 0.5*(3.462435334542 + 3.798170069313) + 7.3, -(0.5*(-6.14367530233 -7.067145485502) - 2.0), 0.0 });

    auto localCorners = corners;
    for (int i = 4; i < 8; ++i)
        localCorners[i][2] += h;

    const auto ab = localCorners[5] - localCorners[4];
    const auto ac = localCorners[6] - localCorners[4];
    auto normal = Dumux::crossProduct(ab, ac);
    normal /= normal.two_norm();

    const auto dist0 = std::max(0.0, (localCorners[4] - ref0)*normal);
    const auto dist1 = std::max(0.0, (localCorners[4] - ref1)*normal);

    const bool dry0 = dist0 < 0.75;
    const bool dry1 = dist1 < 0.75;

    std::cout << "Water heights:      " << dist0 << ", " << dist1 << " --> dry?: " << dry0 << ", " << dry1 << std::endl;

    const auto rotateX = make3DRotation({1.0, 0.0, 0.0}, angle[0]);
    const auto rotateY = make3DRotation({0.0, 1.0, 0.0}, angle[1]);
    const auto elevation0 = rotateX(rotateY(ref0))[2];
    const auto elevation1 = rotateX(rotateY(ref1))[2];

    std::cout << "Channel elevations: " << elevation0 << ", " << elevation1 << std::endl;

    const auto p0 = 1000*9.81*(dist0 + elevation0)*1e-3;
    const auto p1 = 1000*9.81*(dist1 + elevation1)*1e-3;
    std::cout << "Total pressure:     " << p0 << ", " << p1 << std::endl;
    std::cout << "========================================================================" << std::endl;

    {
        std::ofstream metaData(outputName + ".txt");
        metaData << Dumux::Fmt::format("{} {} {} {} {} {} {}\n", angle[0], angle[1], volume, p0, p1, dry0, dry1);
    }

    return std::array<std::pair<double, bool>, 2>({
        std::make_pair( p0, dry0 ), std::make_pair( p1, dry1 )
    });
}

template<class Point, class Tree>
std::array<double, 2> pressureGradients(const Tree& aabbTree, const std::vector<Point>& corners, std::array<double, 2> angle, std::array<double, 2> volumes, int timeStepIndex)
{
    // for the first reservoir use the given angles
    const auto state0 = computeChannelStates(aabbTree, corners, angle, volumes[0], timeStepIndex, 0);
    // for the second reservoir we use minus the given angles (geometry mirrored at center of rotation (origin))
    const auto state1 = computeChannelStates(aabbTree, corners, { -angle[0], -angle[1] }, volumes[1], timeStepIndex, 1);
    // the pressures are also mirrorred so they have to be compared with the respective opposite pressure (TODO check sign?)
    auto diffP0 = state0[0].first - state1[1].first;
    // flow is only happening if the upstream connection is not dry
    if ((diffP0 > 0.0 && state0[0].second) || (diffP0 < 0.0 && state1[1].second))
        diffP0 = 0.0;
    auto diffP1 = state0[1].first - state1[0].first;
    // flow is only happening if the upstream connection is not dry
    if ((diffP1 > 0.0 && state0[1].second) || (diffP1 < 0.0 && state1[0].second))
        diffP1 = 0.0;

    if (diffP0 > 1e-7)
        std::cout << "Channel (1) flow direction R0 " << std::string(int(std::abs(diffP0)), '>') << "> R1 (ΔP: " << diffP0 << " Pa)\n";
    else if (diffP0 < -1e-7)
        std::cout << "Channel (1) flow direction R0 <" << std::string(int(std::abs(diffP0)), '<') << " R1 (ΔP: " << diffP0 << " Pa)\n";
    else
        std::cout << "Channel (1) flow stopped   R0 |-----| R1\n";

    if (diffP1 > 1e-7)
        std::cout << "Channel (2) flow direction R0 " << std::string(int(std::abs(diffP1)), '>') << "> R1 (ΔP: " << diffP1 << " Pa)\n";
    else if (diffP1 < -1e-7)
        std::cout << "Channel (2) flow direction R0 <" << std::string(int(std::abs(diffP1)), '<') << " R1 (ΔP: " << diffP1 << " Pa)\n";
    else
        std::cout << "Channel (2) flow stopped   R0 |-----| R1\n";

    return { diffP0, diffP1 };
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    Parameters::init(argc, argv);

    using ALU = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;
    Dumux::GridManager<ALU> gridManagerAlu;
    gridManagerAlu.init();
    const auto gridView = gridManagerAlu.grid().leafGridView();

    Dune::VTKWriter<ALU::LeafGridView> vtkWriter(gridView);
    vtkWriter.write("grid", Dune::VTK::base64);

    using EntitySet = GridViewGeometricEntitySet<ALU::LeafGridView, 0>;
    auto aabbTree = std::make_shared<BoundingBoxTree<EntitySet>>(std::make_shared<EntitySet>(gridView));

    double volume = 0.0;
    using Point = Dune::FieldVector<double, 3>;
    Point lowerLeft(1e100), upperRight(-1e100);
    for (const auto& element : elements(gridView))
    {
        const auto geometry = element.geometry();
        volume += geometry.volume();

        for (int i = 0; i < geometry.corners(); ++i)
        {
            const auto& corner = geometry.corner(i);
            for (int dir = 0; dir < 3; ++dir)
            {
                lowerLeft[dir] = std::min(lowerLeft[dir], corner[dir]);
                upperRight[dir] = std::max(upperRight[dir], corner[dir]);
            }
        }
    }

    const auto height = upperRight[2]-lowerLeft[2];
    lowerLeft[2] -= height;
    upperRight[2] += height;

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

    // using HexGeometry = Dune::MultiLinearGeometry<double, 3, 3>;
    // const auto hexCut = HexGeometry(Dune::GeometryTypes::cube(3), corners);
    // writeHex(hexCut, "hex");

    // Dune::Timer timer;
    // const auto treeIntersections = intersectingEntities(hexCut, *aabbTree);
    // std::cout << "Computed " << treeIntersections.size() << " tree intersections in " << timer.elapsed() << std::endl;

    // const auto intersections = convertIntersections<Point>(treeIntersections);
    // writeIntersections<Point>(intersections, "intersections");

    // const double intersectionVolume = computeVolume<Point>(intersections);

    // std::cout << "Total volume: " << volume << " μl" << std::endl;
    // std::cout << "Intersection volume: " << intersectionVolume << " μl (" << intersectionVolume/volume*100 << "%)" << std::endl;

    // rotate once around the clock at a given angle and compute the water body geometry for constant volume
    const auto angleDegree = getParam<double>("Problem.Angle", 18.0);
    const auto theta = angleDegree/180.0*M_PI;
    const auto sinTheta = std::sin(theta);

    const auto rotationsPerSecond = getParam<double>("Problem.RotationsPerMinute", 4.0)/60.0;
    const auto cycles = getParam<double>("TimeLoop.Cycles", 1.0);
    const auto tEnd = getParam<double>("TimeLoop.TEnd", cycles/rotationsPerSecond);
    const auto dt = getParam<double>("TimeLoop.Dt", tEnd/cycles/200.0);
    auto timeLoop = std::make_shared<TimeLoop<double>>(0.0, dt, tEnd);

    const auto initialVolume = getParam<double>("Problem.InitialVolumeInMicroLiter", 300.0);
    const auto channelVolume = getParam<double>("Problem.SingleChannelVolumeInMicroLiter", 16.735);
    std::array<double, 2> volumes({ 0.0, initialVolume - 2*channelVolume });

    // m^3/(s*Pa)
    const auto channelTransmissibility = getParam<double>("Problem.ChannelTransmissibility", 9e-10);
    // m/(s*Pa)
    const auto channelVelFactor = getParam<double>("Problem.ChannelVelFactor", 0.026e-1);

    const auto outputFileName = getParam<std::string>("Problem.OutputFileName", "output.txt");
    std::ofstream output(outputFileName);
    output << "Time[s] volTotal[μl] volA[μl] volB[μl] flux_ch0[μl/s] flux_ch1[μl/s] maxv_ch0[m/s] maxv_ch1[m/s] beta[rad] gamma[rad]\n";
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
        // regularized volume so that the bracket algorithm for finding water heights is guaranteed to work
        auto volRegu = volumes;
        volRegu[0] = std::max(1e-5, volRegu[0]);
        volRegu[1] = std::max(1e-5, volRegu[1]);
        const auto diffP = pressureGradients(*aabbTree, corners, { gamma, beta }, volRegu, timeLoop->timeStepIndex());

        // Update volumes using pressure gradients and channel resistance
        // convert to μl/s
        const auto netFluxPredict = -1e9*(diffP[0]*channelTransmissibility + diffP[1]*channelTransmissibility);
        // make sure only as much flows as is actually there
        const auto netVolumeExchange = std::clamp(netFluxPredict*timeLoop->timeStepSize(),
            std::max(-volumes[0], volumes[1]-volume), // minimum is contrained by what's left in volume 0 or can be added to volume 1
            std::min(volumes[1], volume-volumes[0]) // minimum is contrained by what's left in volume 1 or can be added to volume 0
        );

        volumes[0] += netVolumeExchange;
        volumes[1] -= netVolumeExchange;

        std::cout << "Time: " << curTime << ", t: " << t << std::endl;
        output << curTime << " "
               << volumes[0] + volumes[1] + 2*channelVolume << " "
               << volumes[0] << " "
               << volumes[1] << " "
               << -1e9*(diffP[0]*channelTransmissibility) << " "
               << -1e9*(diffP[1]*channelTransmissibility) << " "
               << -(diffP[0]*channelVelFactor) << " "
               << -(diffP[1]*channelVelFactor) << " "
               << beta << " "
               << gamma
               << "\n";

        // go to next time step
        timeLoop->advanceTimeStep();

        // report statistics of this time step
        timeLoop->reportTimeStep();
        std::cout << std::endl;

    } while (!timeLoop->finished());

    timeLoop->finalize(gridView.comm());

    return 0;
}