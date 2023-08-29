#include <iostream>
#include "integration.hpp"
// #include "tetrahedron_jaskowiec_rule.hpp"
#include "mesh.hpp"
#include <set>
#include <chrono>

using namespace GaussLobatto;
using namespace IntegrationMonomial;
using namespace Gauss;

int main()
{
    /*
        // Test [2-9]-points Gauss-Lobatto rule in [-1,1]
        for (int n = 2; n < 10; n++)
        {
            auto p = computeGaussLobattoXW(n);
            for (const auto &X : p.first)
            {
                std::cout << X << " ";
            }
            std::cout << std::endl;
            real sum = 0.0;
            for (const auto &W : p.second)
            {
                std::cout << W << " ";
                sum += W;
            }
            std::cout << std::endl
                      << sum << std::endl;
        }

        Point3D point0(0, 0, 0);
        Point3D point1(1, 0, 0);
        auto edge = Edge(point0, point1);
        std::cout << edge << edge[0] << edge[1] << std::endl;

        // Test n-points Gauss-Lobatto points and weights for an edge
        auto p = computeGaussLobattoPointsWOnEdge(edge, 3);
        for (const auto &point : p.first)
        {
            std::cout << point << " ";
        }
        std::cout << std::endl;
        for (const auto &weight : p.second)
        {
            std::cout << weight << " ";
        }
        std::cout << std::endl;
    */

    // Test quadrature-free integration of a monomial
    std::string filename = "test.geo";
    Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> mesh(filename);

    for (std::size_t p = 0; p < 1; p++)
    {

        Polyhedron<Polygon3D> P = mesh.getPolyhedron(p);
        /*
                std::set<Point3D> vertices;
                for (std::size_t f = 0; f < P.numPolygons(); f++)
                {
                    for (std::size_t e = 0; e < P[f].numEdges(); e++)
                    {
                        // std::cout<<p1[f][e][0]<<" "<<p1[f][e][1]<<std::endl;
                        vertices.insert(P[f][e][0]);
                    }
                }
                // Print vertices
                for (const auto &v : vertices)
                {
                    std::cout << v << std::endl;
                }
                std::cout << std::endl;
                */
        /*
                // Test integration over polyhedra of monomials embedded in space, global reference system
                // Compute integrals of all monomials in 3D up to order N
                std::size_t N = 4;
                std::vector<Monomial3D> monomialOrdered3D = Monomial3D::getMonomialsOrdered(N);
                for (std::size_t i = 0; i < monomialOrdered3D.size(); i++)
                {
                    // Measure the execution time of the integrateMonomial function
                    auto start = std::chrono::high_resolution_clock::now();
                    real result = integrateMonomial(3, P, monomialOrdered3D[i]);
                    auto end = std::chrono::high_resolution_clock::now();

                    // Calculate the elapsed time
                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

                    // Print the result and execution time
                    std::cout << "Integral of " << monomialOrdered3D[i] << ": " << result << std::endl;
                    std::cout << "Execution time: " << duration << " microseconds" << std::endl
                              << std::endl;
                }
        */

        // Test integration over polyhedra of monomials embedded in space, local reference system
        // Compute integrals of all monomials in 3D up to order N
        std::size_t N = 1;
        std::vector<Monomial3D> monomialOrdered3D = Monomial3D::getMonomialsOrdered(N);
        Point3D X_P = getPolyhedronCentroid(P);
        // Point3D X_P(1.01786, 1.35714, 0.75);
        auto h_P = P.getDiameter();
        for (std::size_t i = 0; i < monomialOrdered3D.size(); i++)
        {
            // Measure the execution time of the integrateMonomial function
            auto start = std::chrono::high_resolution_clock::now();
            real result = integrateMonomial(3, P, monomialOrdered3D[i], X_P, h_P);
            auto end = std::chrono::high_resolution_clock::now();

            // Calculate the elapsed time
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            // Print the result and execution time
            std::cout << "Integral of " << monomialOrdered3D[i] << ": " << result << std::endl;
            std::cout << "Execution time: " << duration << " microseconds" << std::endl
                      << std::endl;
        }
    }
    std::cout << std::endl
              << std::endl;
    /*
        Polygon3D F = mesh.getPolygon(1);

        std::set<Point3D> vertices;
        for (std::size_t e = 0; e < F.numEdges(); e++)
        {
            // std::cout<<p1[f][e][0]<<" "<<p1[f][e][1]<<std::endl;
            vertices.insert(F[e][0]);
        }
        // Print vertices
        for (const auto &v : vertices)
        {
            std::cout << v << std::endl;
        }
        std::cout << std::endl;
    */
    Point3D p0(0, 0, 0);
    Point3D p1(3, 0, 0);
    Point3D p2(3, 2, 0);
    Point3D p3(1.5, 4, 0);
    Point3D p4(0, 4, 0);

    Edge3D e1(p0, p1);
    Edge3D e2(p1, p2);
    Edge3D e3(p2, p3);
    Edge3D e4(p3, p4);
    Edge3D e5(p4, p0);

    Point3D p10(2, 0, 0);
    Point3D p20(0, 2, 0);

    Edge3D e10(p0, p10);
    Edge3D e20(p10, p20);
    Edge3D e30(p20, p0);

    // Polygon3D F{e1, e2, e3, e4, e5};
    auto F = mesh.getPolygon(1);

    Point3D X_F = getPolygonCentroid(F);
    real h_F = F.getDiameter();
    Point3D e_x = F.get_e_x();
    Point3D e_y = F.get_e_y();
    real A_F = F.getArea();
    std::cout << "Centroid = " << X_F << std::endl
              << "Diameter = " << h_F << std::endl
              << "First local axis = " << e_x << std::endl
              << "Second local axis = " << e_y << std::endl
              << "Area = " << A_F << std::endl
              << std::endl;

    std::cout << integrateMonomial(2, F, Monomial2D(2, 0, 1.0), X_F, h_F, e_x, e_y) / A_F << std::endl;

    // std::cout << integrateMonomial(2, F, Monomial2D(2, 0, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    /*
    std::cout << integrateMonomial(2, F, Monomial2D(3, 0, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(4, 0, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(5, 0, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(6, 0, 1.0), X_F, h_F, e_x, e_y) << std::endl
              << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(0, 1, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(0, 2, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(0, 3, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(0, 4, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(0, 5, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(0, 6, 1.0), X_F, h_F, e_x, e_y) << std::endl
              << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(1, 1, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(2, 1, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(3, 1, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(3, 2, 1.0), X_F, h_F, e_x, e_y) << std::endl
              << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(1, 1, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(1, 2, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(1, 3, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    std::cout << integrateMonomial(2, F, Monomial2D(2, 3, 1.0), X_F, h_F, e_x, e_y) << std::endl;
*/

    // Test integration over polygons of monomials embedded in plane, global reference system
    // Compute integrals of all monomials in 2D up to order N
    std::size_t N = 20;
    std::vector<Monomial2D> monomialOrdered2D = Monomial2D::getMonomialsOrdered(N);
    /*
    for (std::size_t i = 0; i < monomialOrdered2D.size(); i++)
    {
        // Measure the execution time of the integrateMonomial function
        auto start = std::chrono::high_resolution_clock::now();
        real result = integrateMonomial(2, F, monomialOrdered2D[i]);
        auto end = std::chrono::high_resolution_clock::now();

        // Calculate the elapsed time
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        // Print the result and execution time
        std::cout << "Integral of " << monomialOrdered2D[i] << ": " << result << std::endl;
        std::cout << "Execution time: " << duration << " microseconds" << std::endl
                  << std::endl;
    }

    // Test integration over polygons of monomials embedded in plane, local reference system
    //MonomialsFaceIntegralsCache::initialize(F,20);
    for (std::size_t i = 0; i < monomialOrdered2D.size(); i++)
    {
        // Measure the execution time of the integrateMonomial function
        auto start1 = std::chrono::high_resolution_clock::now();
        real result1 = integrateMonomial(2, F, monomialOrdered2D[i], X_F, h_F, e_x, e_y);
        auto end1 = std::chrono::high_resolution_clock::now();

        auto start2 = std::chrono::high_resolution_clock::now();
        real result2 = MonomialsFaceIntegralsCache::getCache(F,monomialOrdered2D[i]);
        auto end2 = std::chrono::high_resolution_clock::now();

        // Calculate the elapsed time
        auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1).count();
        auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2).count();

        // Print the result and execution time
        std::cout << "Integral of " << monomialOrdered2D[i] << ": " << result1 << std::endl;
        std::cout << "Execution time: " << duration1 << " microseconds" << std::endl;
        std::cout << "Integral of " << monomialOrdered2D[i] << ": " << result2 << std::endl;
        std::cout << "Execution time: " << duration2 << " microseconds" << std::endl
                  << std::endl;
    }
*/
    /*
     // Test integration over polygons of monomials embedded in space restricted to polygons
     std::cout << "Negative z face (bottom)" << std::endl;
     std::cout << integrateMonomial(2, P[0], Monomial3D::getMonomialsOrdered(2)[0]) << std::endl;
     std::cout << integrateMonomial(2, P[0], Monomial3D::getMonomialsOrdered(2)[1]) << std::endl;
     std::cout << integrateMonomial(2, P[0], Monomial3D::getMonomialsOrdered(2)[4]) << std::endl;
     std::cout << integrateMonomial(2, P[0], Monomial3D::getMonomialsOrdered(2)[5]) << std::endl;
     std::cout << "Positive x face" << std::endl;
     std::cout << integrateMonomial(2, P[1], Monomial3D::getMonomialsOrdered(2)[0]) << std::endl;
     std::cout << integrateMonomial(2, P[1], Monomial3D::getMonomialsOrdered(2)[1]) << std::endl;
     std::cout << integrateMonomial(2, P[1], Monomial3D::getMonomialsOrdered(2)[4]) << std::endl;
     std::cout << integrateMonomial(2, P[1], Monomial3D::getMonomialsOrdered(2)[5]) << std::endl;
     */
    /*
        //MonomialsFaceIntegralsCache::initialize(mesh,20);
        // Test to integrate monomials embedded in space restricted to polygons
        // Compute integrals of all monomials in 3D up to order 4
        std::vector<Monomial3D> monomialOrdered3D = Monomial3D::getMonomialsOrdered(10);
        Point3D X_P(0.5, 0.5, 0.5);
        real h_P = sqrt(3);
        for (std::size_t i = 0; i < monomialOrdered3D.size(); i++)
        {
            // Measure the execution time of the integrateMonomial function
            auto start = std::chrono::high_resolution_clock::now();
            real result = integrateMonomial3DRestrictedMonomial2D(X_P, h_P, F, monomialOrdered3D[i], monomialOrdered2D[0]);
            auto end = std::chrono::high_resolution_clock::now();

            // Calculate the elapsed time
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            // Print the result and execution time
            std::cout << "Integral on the face F of " << monomialOrdered3D[i] << ": " << result << std::endl;
            std::cout << "Execution time: " << duration << " microseconds" << std::endl
                      << std::endl;
        }
    */

    // Test integrateTetrahedron
    Point3D v1(0.1, 5, 2);
    Point3D v2(1. / 3., 1.2, 31);
    Point3D v3(0.1, 4. / 5., 0.01);
    Point3D v4(5, 6.7, 92);

    // auto XYZW = gaussPointsTetrahedron(v1, v2, v3, v4, 1);

    auto func1 = [](real x, real y, real z)
    {
        return x * y;
    };
    auto I = integrateFunctionOverTetrahedron(func1, v1, v2, v3, v4, 2);
    std::cout << I << std::endl;

    // Test integratePolyhedron
    auto P = mesh.getPolyhedron(0);
    Point3D X_P = getPolyhedronCentroid(P);
    real h_P = P.getDiameter();
    auto m = Monomial3D::getMonomialsOrdered(3)[4];
    auto func0 = [&X_P, &h_P, &m](real x, real y, real z)
    {
        return std::pow((x - X_P[0]) / h_P, 7);
    };
    auto func = [&X_P, &h_P, &func0, &m](real x, real y, real z)
    {
        return func0(x,y,z) * m.evaluate((Point3D(x, y, z) - X_P) / h_P);
        // return x;
    };
    std::cout<<std::setprecision(16);
    auto start = std::chrono::high_resolution_clock::now();
    I = integrateFunctionOverPolyhedron(func, P, 5);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_t = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Integral result = " << I << ", elapsed time = " << elapsed_t << std::endl;
    start = std::chrono::high_resolution_clock::now();
    I = integrateMonomial(3, P, Monomial3D(9, 0, 0, 1.0), X_P, h_P);
    end = std::chrono::high_resolution_clock::now();
    elapsed_t = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Integral result = " << I << ", elapsed time = " << elapsed_t << std::endl;

    return 0;
}
/*
for(std::size_t n=4;n<9;n++)
{
    std::size_t size = (n - 3) / 2 +1;
    std::cout<<n<<" "<<size<<std::endl;
}
*/
/*
#include "integration.hpp"
#include <iostream>

int main()
{
    // create three points
    Point3D point1(-1, 0, 0);
    Point3D point2(1, 0, 0);
    Edge edge1(point1, point2);

    int n = 5; // Order of the Legendre polynomial
    int numRoots = n;
    double tol = 1e-9;

    for (int i = 1; i <= numRoots; ++i)
    {
        // Create a functor representing the derivative of the Legendre polynomial of order n
        auto func = [n](double x)
        { return legendreDerivative(n, x); };

        // Find the root using toms748_solve
        boost::uintmax_t maxIterations = 100;
        auto result = boost::math::tools::toms748_solve(func, -1.0, 1.0, tol, maxIterations);

        double root = result.first; // Extract the root from the result pair
        double estimatedError = result.second; // Extract the estimated error from the result pair
        std::cout<<root<<" ";
    }

    /*
        for (int order = 2; order < 3; order++)
        {
            auto gausspoints = computeGaussLobattoPointsOnEdge(edge1, order);
            for (const auto &gausspoint : gausspoints)
            {
                std::cout << gausspoint << " ";
            }
            std::cout<<std::endl;
        }

    return 0;
}
*/