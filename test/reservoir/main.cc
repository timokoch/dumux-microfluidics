
#include <config.h>

#include <vector>

#include <dune/common/timer.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/common/math.hh>
#include <dumux/nonlinear/findscalarroot.hh>

#include <dumux/io/grid/gridmanager_alu.hh>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>

#include <test/geometry/writetriangulation.hh>

template<class Point, class TreeIntersections>
auto convertIntersections(const TreeIntersections& treeIntersections)
{
    Dune::Timer timer;
    std::vector<std::vector<Point>> intersections;
    intersections.reserve(treeIntersections.size());
    for (const auto& is : treeIntersections)
        intersections.emplace_back(std::vector<Point>(is.corners()));
    std::cout << "Converted to output format in " << timer.elapsed() << " seconds." << std::endl;
    return intersections;
}

template<class Point>
void writeIntersections(const std::vector<std::vector<Point>>& intersections, const std::string& name)
{
    Dune::Timer timer;
    std::cout << "Writing " << intersections.size() << " intersections to file ...";
    Dumux::writeVTUTetrahedron(intersections, name);
    std::cout << " done ( " << timer.elapsed() << " seconds)." << std::endl;
}

template<class Point>
void writeIntersections(const std::vector<std::array<Point, 4>>& intersections, const std::string& name)
{
    Dune::Timer timer;
    std::cout << "Writing " << intersections.size() << " intersections to file ...";
    Dumux::writeVTUTetrahedron(intersections, name);
    std::cout << " done ( " << timer.elapsed() << " seconds)." << std::endl;
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

template<class Point, class Tree>
void optimizeVolume(const Tree& aabbTree, std::vector<Point> corners, std::array<double, 2> angle, double volume, int index)
{
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
    auto treeIntersections = intersectingEntities(hexCut, aabbTree);
    auto intersections = convertIntersections<Point>(treeIntersections);

    // int count = 0;
    const auto residual = [&](const double h)
    {
        auto localCorners = corners;
        for (int i = 4; i < 8; ++i)
            localCorners[i][2] += h;

        hexCut = HexGeometry(Dune::GeometryTypes::cube(3), localCorners);
        treeIntersections = intersectingEntities(hexCut, aabbTree);
        intersections = convertIntersections<Point>(treeIntersections);
        const auto iVol = computeVolume<Point>(intersections);
        // writeIntersections<Point>(intersections, "debug-intersections-" + std::to_string(index) + "_" + std::to_string(count));
        // writeHex(hexCut, "debug-hex-" + std::to_string(index) + "_" + std::to_string(count++));

        // if (index == 10 && count == 9)
        // {
        //     using IntersectionAlgorithm = Dumux::GeometryIntersection<HexGeometry, HexGeometry>;
        //     using Intersection = typename IntersectionAlgorithm::Intersection;
        //     Intersection intersection;

        //     if (IntersectionAlgorithm::intersection(hexCut, tet, intersection, true))
        //     {
        //         const auto triangulation = Dumux::triangulate<3, 3>(intersection);
        //         writeIntersections<Point>(triangulation, "debug-debug-intersections-" + std::to_string(index) + "_" + std::to_string(count));
        //     }
        //     else
        //     {
        //         std::cout << "No intersection found!!" << std::endl;
        //     }
        // }

        std::cout << "h: " << h << "-> vol: " << iVol << " / residual: " << iVol - volume << std::endl;
        return iVol - volume;
    };

    auto h = Dumux::findScalarRootBrent(-1.0, 1.0, residual, 1e-2);
    auto vol = residual(h) + volume;

    std::cout << "Optimized volume at " << vol << " mm^3" << std::endl;
    writeIntersections<Point>(intersections, "intersections-" + std::to_string(index));

    Point ref0({ 0.5*(3.462435334542 + 3.798170069313) + 7.3, 0.5*(-6.14367530233 -7.067145485502) - 2.0, 0.0 });
    Point ref1({ 0.5*(3.462435334542 + 3.798170069313) + 7.3, -(0.5*(-6.14367530233 -7.067145485502) - 2.0), 0.0 });

    const auto ab = corners[5] - corners[4];
    const auto ac = corners[6] - corners[4];
    auto normal = Dumux::crossProduct(ab, ac);
    normal /= normal.two_norm();

    const auto dist0 = (corners[4] - ref0)*normal;
    const auto dist1 = (corners[4] - ref1)*normal;

    std::cout << "Water heights: " << dist0 << ", " << dist1 << std::endl;
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
        corners[i][2] -= 1.8*height;

    using HexGeometry = Dune::MultiLinearGeometry<double, 3, 3>;
    const auto hexCut = HexGeometry(Dune::GeometryTypes::cube(3), corners);
    writeHex(hexCut, "hex");

    Dune::Timer timer;
    const auto treeIntersections = intersectingEntities(hexCut, *aabbTree);
    std::cout << "Computed " << treeIntersections.size() << " tree intersections in " << timer.elapsed() << std::endl;

    const auto intersections = convertIntersections<Point>(treeIntersections);
    writeIntersections<Point>(intersections, "intersections");

    const double intersectionVolume = computeVolume<Point>(intersections);

    std::cout << "Total volume: " << volume << " mm^3" << std::endl;
    std::cout << "Intersection volume: " << intersectionVolume << " mm^3 (" << intersectionVolume/volume*100 << "%)" << std::endl;

    // rotate once around the clock at a given angle and compute the water body geometry for constant volume
    const auto angleDegree = getParam<double>("Angle", 18.0);
    const std::size_t samples = getParam<std::size_t>("Samples", 36);
    const auto theta = angleDegree/180.0*M_PI;
    const auto sinTheta = std::sin(theta);
    const auto t = Dumux::linspace(0.0, 2*M_PI, samples);
    for (int i = 0; i < t.size(); ++i)
    {
        const auto sinT = std::sin(t[i]);
        const auto cosT = std::cos(t[i]);
        const auto gamma = std::asin(-sinTheta*sinT);
        const auto beta = std::asin(-sinTheta*cosT/std::sqrt(1.0 - sinTheta*sinTheta*sinT*sinT));

        optimizeVolume(*aabbTree, corners, { gamma, beta }, intersectionVolume, i);
    }

    return 0;
}