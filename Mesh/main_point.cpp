#include "point.hpp"
#include "traits.hpp"
#include <iostream>

using namespace geometry;

int main()
{
    // Create 3D points
    Point<real, real, real> p1(1, 0.0, 0);
    Point<real, real, real> p2(0, 1, 0);

    std::cout << p1[0] << std::endl;
    std::cout << "dimension = " << p1.getDimension() << std::endl;
    std::cout << "coordinates =" << p1.getCoordinates()[0] << std::endl;

    Point p3(p1);
    // p3.setCoordinates(0, 1, 0);
    Point<real, real> p2D(1, 2);

    // Compute the sum of two points using the + operator
    auto sumPoint = p1 + p2;
    std::cout << "Sum Point: (" << sumPoint[0] << ", "
              << sumPoint[1] << ", " << sumPoint[2] << ")" << std::endl;

    // Compute the difference of two points using the - operator
    auto diffPoint = p1 - p2;
    std::cout << "Difference Point: (" << diffPoint[0] << ", "
              << diffPoint[1] << ", " << diffPoint[2] << ")" << std::endl;

    // Compute the piecewise multiplication between two points using inclass method
    auto multPoint = p1.piecewiseMultiply(p2);
    std::cout << "piecewise mult: (" << multPoint[1] << ", "
              << multPoint[2] << ")" << std::endl;

    /*
        // Compute the piecewise multiplication between two points using inclass method
        auto multPoint1 = piecewiseMultiply(p1,p2);
        std::cout << "piecewise mult: (" << multPoint1[1] << ", "
                  << multPoint1[2] << ")" << std::endl;
    */

    // distance
    auto dist_inClass = p1.distance(p2);
    auto dist_outClass = distance(p1, p2);
    std::cout << "distance in class: " << dist_inClass << std::endl;
    std::cout << "distance out of class: " << dist_outClass << std::endl;

    // scalar multiplication
    auto scalarTimesPoint = 2 * p3;
    auto PointTimesScalar = p3 * 2;
    std::cout << "scalarTimesPoint: (" << scalarTimesPoint[0] << ", "
              << scalarTimesPoint[1] << ", " << scalarTimesPoint[2] << ")" << std::endl;
    std::cout << "PointTimesScalar: (" << PointTimesScalar[0] << ", "
              << PointTimesScalar[1] << ", " << PointTimesScalar[2] << ")" << std::endl;

    // dot product
    auto dotprod_inClass = p2.dot(p3);
    auto dotprod_outClass = dot(p2, p3);
    std::cout << "dot product in class: " << dotprod_inClass << std::endl;
    std::cout << "dot product out of class: " << dotprod_outClass << std::endl;

    // cross product
    auto crossprod_inClass = p2.cross(p3);
    auto crossprod_outClass = cross(p2, p3);
    std::cout << "cross product in class: "
              << "(" << crossprod_inClass[0] << ", "
              << crossprod_inClass[1] << ", "
              << crossprod_inClass[2] << ")" << std::endl;
    std::cout << "cross product in class: "
              << "(" << crossprod_outClass[0] << ", "
              << crossprod_outClass[1] << ", "
              << crossprod_outClass[2] << ")" << std::endl;

    // compare two points
    if (p1 == p2)
    {
        std::cout << "p1=p2" << std::endl;
    }
    else
    {
        std::cout << "p1!=p2" << std::endl;
    }

    return 0;
}