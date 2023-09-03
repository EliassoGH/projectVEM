#include "traits.hpp"
#include "mesh.hpp"
#include <iostream>
#include <string>
#include <vector>

using namespace geometry;

int main()
{
    std::string filename = "N2.geo";
    Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> mesh(filename);

    std::cout << mesh.numVertices() << std::endl;
    for (std::size_t v = 0; v < mesh.numVertices(); v++)
    {
        std::cout << mesh.getVertex(v) << std::endl;
    }

    std::cout << mesh.numEdges() << std::endl;
    for (std::size_t e = 1; e <= mesh.numEdges(); e++)
    {
        std::cout << mesh.getEdge(e) << " " << mesh.getEdge(-e) << std::endl;
    }

    std::cout << mesh.numPolygons() << std::endl;
    for (std::size_t f = 1; f <= mesh.numPolygons(); f++)
    {
        std::cout << mesh.getPolygon(f) << " " << mesh.getPolygon(-f) << std::endl;
    }

    std::cout << mesh.numPolyhedra() << std::endl;
    for (std::size_t p = 0; p < mesh.numPolyhedra(); p++)
    {
        std::cout << mesh.getPolyhedron(p) << std::endl;
    }

    mesh.print();

    return 0;
}