#include "traits.hpp"
#include "point.hpp"
#include "edge.hpp"
#include <iostream>

using namespace geometry;

int main()
{
    // create three points
    Point3D point1(0, 0, 0);
    Point3D point2(1, 0, 0);
    Point3D point3(0, 1, 0);
    std::cout<<point1<<", "<<point2<<", "<<point3<<std::endl;

    // test copy constructor
    Point3D point4(point2);
    std::cout<<point4<<std::endl;

    // create edge, no need to specify Edge<Point3D> since point2 and point3 are Point3D
    // (C++17 class template argument deduction)
    Edge edge1(point1, point2);
    std::cout << "first point of edge1 is " << edge1[0] << ", second point of edge1 is "
              << edge1[1] << std::endl;

    // test assignment
    //point1 = point3;
    std::cout<<point1<<std::endl;
    edge1.update();
    std::cout << "after assignment first point of edge1 is " << edge1[0] << ", second point of edge1 is "
              << edge1[1] << std::endl;
    std::cout << "now check what happened to the other point which was part of edge1 before: "
              << point1 << std::endl
              << std::endl;

    // test create twin halfedge with same points
    Edge edge2(point1, point2, true);
    std::cout << "first point is " << edge2[0] << ", second point is "
              << edge2[1] << std::endl;

    // equivalence
    std::cout << "the 2 edges are equal? " << (edge1 == edge2) << std::endl;
    std::cout << "edge1 < edge2? " << (edge1 < edge2) << std::endl;

    // properties
    std::cout << "length is " << edge1.getLength() << std::endl;
    std::cout << "direction is " << edge1.getDirection() << std::endl;
    std::cout << "length of twin is " << edge2.getLength() << std::endl;
    std::cout << "direction of twin is " << edge2.getDirection() << std::endl;

    // test createPairEdges
    auto edgesPair=createPairEdges(point3, point2,10);
    std::cout<<edgesPair.first<<", "<<edgesPair.second<<std::endl;
    std::cout<<edgesPair.first.getId()<<", "<<edgesPair.second.getId()<<std::endl;
    std::cout<<edgesPair.first.getOtherHalfEdge()<<", "<<edgesPair.second.getOtherHalfEdge()<<std::endl;

    return 0;
}