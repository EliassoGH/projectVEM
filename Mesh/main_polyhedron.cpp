#include "traits.hpp"
#include "point.hpp"
#include "edge.hpp"
#include "polygon.hpp"
#include "polyhedron.hpp"
#include <iostream>

using namespace geometry;

int main()
{
    // create four points
    Point3D point1(0, 0, 0);
    Point3D point2(1, 0, 0);
    Point3D point3(1, 1, 0);
    Point3D point4(0, 1, 0);
    Point3D point5(0, 0, 1);
    Point3D point6(1, 0, 1);
    Point3D point7(1, 1, 1);
    Point3D point8(0, 1, 1);

    // create edges
    Edge3D edge1(point1, point2);
    Edge3D edge2(point2, point3);
    Edge3D edge3(point3, point4);
    Edge3D edge4(point4, point1);

    Edge3D edge5(point2, point6);
    Edge3D edge6(point6, point5);
    Edge3D edge7(point5, point1);

    Edge3D edge8(point6, point3);
    Edge3D edge9(point4, point5);

    Edge3D edge60(point5, point6);

    // create polygons
    Polygon3D polygon1{edge1, edge2, edge3, edge4};
    Polygon3D polygon2{edge1, edge5, edge6, edge7};
    Polygon3D polygon3{};
    polygon3.addEdge(edge60);
    polygon3.addEdge(edge8);
    polygon3.addEdge(edge3);
    polygon3.addEdge(edge9);
    Polygon3D polygon4{edge4, edge3, edge2, edge1}; // equivalent polygon as polygon1, but different id
    std::cout<<"polygon1==polygon4? "<<(polygon1==polygon4)<<std::endl;

    // check consistency
    std::cout << polygon1.areEdgesConsistent() << std::endl;
    std::cout << polygon2.areEdgesConsistent() << std::endl;
    std::cout << polygon3.areEdgesConsistent() << std::endl;

    // get ids
    std::cout << "id of polygon 1 is " << polygon1.getId() << std::endl;
    std::cout << "id of polygon 2 is " << polygon2.getId() << std::endl;
    std::cout << "id of polygon 3 is " << polygon3.getId() << std::endl;

    // create polyhedrons
    Polyhedron<Polygon3D> polyhedron1{polygon1, polygon2, polygon3};
    Polyhedron<Polygon3D> polyhedron2{};
    std::cout << polyhedron1.numPolygons() << std::endl;
    std::cout << polyhedron2.numPolygons() << std::endl;

    // stream operator
    std::cout<<polyhedron1<<std::endl;
    std::cout<<polyhedron2<<std::endl;

    // try to create a polyhedron with two coinciding polygons
    std::cout<<"Now try to create a polyhedron with two coinciding polygons"<<std::endl;
    Polyhedron<Polygon3D> polyhedron3{polygon1, polygon2, polygon4};
    std::cout << polyhedron3.numPolygons() << std::endl;

    return 0;
}