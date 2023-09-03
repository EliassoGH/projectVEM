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
                  << "sum of the weights is " << sum << std::endl;
    }

    // Test quadrature-free integration of a monomial
    std::string filename = "test.geo";
    Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> mesh(filename);

    for (std::size_t p = 0; p < 1; p++)
    {

        Polyhedron<Polygon3D> P = mesh.getPolyhedron(p);

        // Test integration over polyhedra of monomials embedded in space, global reference system
        // Compute integrals of all monomials in 3D up to order N
        std::cout << std::endl
                  << "Test integration over polyhedra of monomials embedded in space, global reference system" << std::endl;
        std::size_t N = 3;
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

        // Test integration over polyhedra of monomials embedded in space, local reference system
        // Compute integrals of all monomials in 3D up to order N
        std::cout << std::endl
                  << "Test integration over polyhedra of monomials embedded in space, local reference system" << std::endl;
        Point3D X_P = getPolyhedronCentroid(P);
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

    // Test integration over polygons of monomials embedded in plane, global reference system
    // Compute integrals of all monomials in 2D up to order N
    std::cout << std::endl
              << "Test integration over polygons of monomials embedded in plane, global reference system" << std::endl;
    std::size_t N = 4;
    std::vector<Monomial2D> monomialOrdered2D = Monomial2D::getMonomialsOrdered(N);

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
    std::cout << std::endl
              << "Test integration over polygons of monomials embedded in plane, local reference system" << std::endl;
    for (std::size_t i = 0; i < monomialOrdered2D.size(); i++)
    {
        // Measure the execution time of the integrateMonomial function
        auto start1 = std::chrono::high_resolution_clock::now();
        real result1 = integrateMonomial(2, F, monomialOrdered2D[i], X_F, h_F, e_x, e_y);
        auto end1 = std::chrono::high_resolution_clock::now();

        auto start2 = std::chrono::high_resolution_clock::now();
        real result2 = MonomialsFaceIntegralsCache::getCacheMonomial(F, monomialOrdered2D[i]);
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

    // Test integrateTetrahedron
    Point3D v1(0.1, 5, 2);
    Point3D v2(1. / 3., 1.2, 31);
    Point3D v3(0.1, 4. / 5., 0.01);
    Point3D v4(5, 6.7, 92);

    auto func1 = [](real x, real y, real z)
    {
        return x * y;
    };
    auto I = integrateFunctionOverTetrahedron(func1, v1, v2, v3, v4, 2);
    std::cout << I << std::endl;

    // Test integratePolyhedron
    auto P = mesh.getPolyhedron(0);
    Point3D X_P = getPolyhedronCentroid(P);
    auto h_P = P.getDiameter();
    auto m = Monomial3D::getMonomialsOrdered(3)[4];
    auto func0 = [&X_P, &h_P, &m](real x, real y, real z)
    {
        return std::pow((x - X_P[0]) / h_P, 7);
    };
    auto func = [&X_P, &h_P, &func0, &m](real x, real y, real z)
    {
        return func0(x, y, z) * m.evaluate((Point3D(x, y, z) - X_P) / h_P);
        // return x;
    };
    std::cout << std::setprecision(16);
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