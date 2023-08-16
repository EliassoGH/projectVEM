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

    std::cout << mesh.getPolyhedra().at(1)[1][2][0] << std::endl;

    for (const auto &edgePair : mesh.getEdges())
    {
        std::cout << edgePair.second << std::endl;
        std::cout << edgePair.second[0] << " " << edgePair.second[1] << std::endl;
    }

    for (const auto &polygonPair : mesh.getPolygons())
    {
        std::cout << polygonPair.second << std::endl;
    }

    std::cout << mesh.getPolyhedron(0).getPolygon(0).getArea() << std::endl;
    std::cout << mesh.getPolygon(1).getOutwardNormal() << std::endl;
    std::cout << mesh.getPolygon(-1).getOutwardNormal() << std::endl;

    std::cout << "positive edge" << std::endl;
    std::cout << mesh.getPolygon(36).getEdge(0) << std::endl;
    std::cout << mesh.getPolygon(36).getEdge(1) << std::endl;
    std::cout << mesh.getPolygon(36).getEdge(2) << std::endl;
    std::cout << mesh.getPolygon(36).getEdge(3) << std::endl;
    std::cout << mesh.getPolygon(36).getPositiveEdge(0) << std::endl;
    std::cout << mesh.getPolygon(36).getPositiveEdge(1) << std::endl;
    std::cout << mesh.getPolygon(36).getPositiveEdge(2) << std::endl;
    std::cout << mesh.getPolygon(36).getPositiveEdge(3) << std::endl;

    // auto myedge=mesh.getPolyhedron(6)[2][3];
    // std::cout<<myedge<<"with points"<<myedge[0]<<myedge[1]<<"and length"<<myedge.getLength()<<std::endl;

    // std::cout<<mesh.getPolygon(6).second<<std::endl;
    // std::cout<<mesh.getPolygon(5).first<<std::endl;

    // auto a=mesh.polyhedra[1][0][0];
    // auto b=mesh.polyhedra[1][2][0];

    // auto a=mesh.getPolyhedron(1)[0][0];
    // auto b=mesh.getPolyhedron(1)[2][0];
    // std::cout << mesh.getPolyhedron(1)[0][0]<<std::endl;
    // std::cout<< mesh.getPolyhedron(1)[2][0] << std::endl;

    /*
    std::cout << mesh.edges.at(1);
    std::cout << mesh.edges.at(1).getOtherHalfEdge();
    std::cout << mesh.polygons.at(1)[0];
    std::cout << mesh.polygons.at(-1)[0];
    std::cout << mesh.polygons.at(3)[2][0];
    std::cout << mesh.polygons.at(2)[1][0];
    std::cout << mesh.polyhedra.at(1);
    std::cout << mesh.polyhedra.at(1)[2];
    std::cout << mesh.polyhedra.at(1)[2][0];
    std::cout << mesh.polyhedra.at(2)[2][0];
    std::cout << mesh.polyhedra.at(3)[3][0];
    */

    return 0;
}