#include "traits.hpp"
#include "point.hpp"
#include "edge.hpp"
#include "polygon.hpp"
#include <iostream>

using namespace geometry;

int main()
{
    // create four points
    Point3D point1(0, 0, 0);
    Point3D point2(1, 0, 0);
    Point3D point3(1, 1, 0);
    Point3D point4(0, 1, 0);

    // create edges
    Edge3D edge1(point1, point2);
    Edge3D edge2(point3, point2, true);
    Edge3D edge3(point3, point4);
    Edge3D edge4(point4, point1);
    Edge3D edge5(point1, point3);
    Edge3D edge6(point2, point1);

    // create polygon
    Polygon3D polygon1{edge1, edge2, edge3, edge4};
    Polygon3D polygon2;
    polygon2.addEdge(edge1);
    polygon2.addEdge(edge2);
    polygon2.addEdge(edge5);
    Polygon3D polygon3{};

    // print points
    for (std::size_t e = 0; e < polygon1.numEdges(); e++)
    {
        std::cout << polygon1[e] << ": " << polygon1[e][0] << ", " << polygon1[e][1] << std::endl;
    }

    // check consistency
    std::cout << polygon1.areEdgesConsistent() << std::endl;
    std::cout << polygon2.areEdgesConsistent() << std::endl;

    // get ids
    std::cout << "id of polygon 1 is " << polygon1.getId() << std::endl;
    std::cout << "id of polygon 2 is " << polygon2.getId() << std::endl;
    std::cout << "id of polygon 3 is " << polygon3.getId() << std::endl;

    // check if it prevents insertion of same edge
    // Polygon3D polygon4{edge1, edge6};

    // test areCyclicPermutations for vectors
    std::vector<Edge3D> v1{edge1, edge2, edge3, edge4};
    std::vector<Edge3D> v2{edge3, edge2, edge1, edge4};
    std::cout << areEquivalentCyclicPermutations(v1, v2) << std::endl;

    // test == operator for polygons
    Polygon3D p1{edge1, edge2, edge3, edge4};
    Polygon3D p2{edge3, edge4, edge1, edge2};
    std::cout << "id of p1 is " << p1.getId() << std::endl;
    std::cout << "id of p2 is " << p2.getId() << std::endl;
    std::cout << "p1==p2? " << (p1 == p2) << std::endl;
    std::cout << "p1<p2? " << (p1 < p2) << std::endl;

    // test output stream operator
    std::cout << p1 << std::endl;

    // Test compute outward normal
    std::cout<<p1.getOutwardNormal()<<std::endl;
    std::cout<<p1.getArea()<<std::endl;

    // Test compute first local axis e_x
    std::cout<<p1.get_e_x()<<std::endl;

    // Test compute second local axis e_y
    std::cout<<p1.get_e_y()<<std::endl;

    // Test to compute the diameter
    std::cout<<p1.getDiameter()<<std::endl;

    // Test to compute the centroid
    //std::cout<<p1.getCentroid()<<std::endl;

    return 0;
}