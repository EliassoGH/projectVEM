#include <iostream>
#include "integration.hpp"
#include "mesh.hpp"
#include <set>
#include <chrono>

using namespace GaussLobatto;
using namespace IntegrationMonomial;

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
    std::string filename = "N2.geo";
    Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> mesh(filename);

    
    for (std::size_t p = 0; p < mesh.numPolyhedra(); p++)
    {
        /*
        Polyhedron<Polygon3D> P = mesh.getPolyhedron(p);

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
        /*
                // Test integration over polyhedra of monomials embedded in space, local reference system
                // Compute integrals of all monomials in 3D up to order N
                std::size_t N = 4;
                std::vector<Monomial3D> monomialOrdered3D = Monomial3D::getMonomialsOrdered(N);
                for (std::size_t i = 0; i < monomialOrdered3D.size(); i++)
                {
                    // Measure the execution time of the integrateMonomial function
                    auto start = std::chrono::high_resolution_clock::now();
                    real result = integrateMonomial(3, P, monomialOrdered3D[i],Point3D(0.5,0.5,0.5),sqrt(3));
                    auto end = std::chrono::high_resolution_clock::now();

                    // Calculate the elapsed time
                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

                    // Print the result and execution time
                    std::cout << "Integral of " << monomialOrdered3D[i] << ": " << result << std::endl;
                    std::cout << "Execution time: " << duration << " microseconds" << std::endl
                              << std::endl;
                }
        */
    }
    std::cout << std::endl
              << std::endl;

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

    Point3D X_F = getPolygonCentroid(F);
    real h_F = F.getDiameter();
    Point3D e_x = F.get_e_x();
    Point3D e_y = F.get_e_y();
    real A_F = F.getArea();
    std::cout << X_F << std::endl
              << h_F << std::endl
              << e_x << std::endl
              << e_y << std::endl
              << A_F << std::endl;

    //std::cout << integrateMonomial(2, F, Monomial2D(1, 0, 1.0), X_F, h_F, e_x, e_y) << std::endl;
    //std::cout << integrateMonomial(2, F, Monomial2D(2, 0, 1.0), X_F, h_F, e_x, e_y) << std::endl;
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
        std::size_t N = 16;
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
        */
        // Test integration over polygons of monomials embedded in plane, local reference system
        for (std::size_t i = 140; i < monomialOrdered2D.size(); i++)
        {
            // Measure the execution time of the integrateMonomial function
            auto start = std::chrono::high_resolution_clock::now();
            real result = integrateMonomial(2, F, monomialOrdered2D[i], X_F, h_F, e_x, e_y);
            auto end = std::chrono::high_resolution_clock::now();

            // Calculate the elapsed time
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            // Print the result and execution time
            std::cout << "Integral of " << monomialOrdered2D[i] << ": " << result << std::endl;
            std::cout << "Execution time: " << duration << " microseconds" << std::endl
                      << std::endl;
        }
    

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
        // Test to integrate monomials embedded in space restricted to polygons
        // Compute integrals of all monomials in 3D up to order 4
        std::vector<Monomial3D> monomialOrdered3D = Monomial3D::getMonomialsOrdered(4);
        Point3D X_P(0.5, 0.5, 0.5);
        real h_P=sqrt(3);
        for (std::size_t i = 0; i < monomialOrdered3D.size(); i++)
        {
            // Measure the execution time of the integrateMonomial function
            auto start = std::chrono::high_resolution_clock::now();
            real result = integrateMonomial3DRestricted(X_P, h_P, F, monomialOrdered3D[i]);
            auto end = std::chrono::high_resolution_clock::now();

            // Calculate the elapsed time
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            // Print the result and execution time
            std::cout << "Integral on the face F of " << monomialOrdered3D[i] << ": " << result << std::endl;
            std::cout << "Execution time: " << duration << " microseconds" << std::endl
                      << std::endl;
        }
    */
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