#include "virtualProjections.hpp"
#include "solver.hpp"
#include <iostream>
#include <math.h>

int main()
{
    std::string filename = "voro8.geo";
    Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> mesh(filename);

    unsigned int k = 2;
    MonomialsFaceIntegralsCache::initialize(mesh, 2 * k - 1);

    real lambda = 1.0;
    real mu = 1.0;
    real E = 2.5;
    real nu = 0.25;
    VirtualDofsCollection DOFS(mesh, k);

    VirtualFaceProjections vf(DOFS, mesh, k);

    LocalVirtualDofsCollection dofs(mesh, DOFS);

    auto funcx = [&lambda, &mu](real x, real y, real z)
    {
        // return -0.1*std::pow(M_PI,2)*((lambda+mu)*cos(M_PI*x)*sin(M_PI*y+M_PI*z)-(lambda+4*mu)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z));
        return 0.2;
    };
    auto funcy = [&lambda, &mu](real x, real y, real z)
    {
        // return -0.1*std::pow(M_PI,2)*((lambda+mu)*cos(M_PI*y)*sin(M_PI*x+M_PI*z)-(lambda+4*mu)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z));
        return 0.0;
    };
    auto funcz = [&lambda, &mu](real x, real y, real z)
    {
        // return -0.1*std::pow(M_PI,2)*((lambda+mu)*cos(M_PI*z)*sin(M_PI*x+M_PI*y)-(lambda+4*mu)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z));
        return 0.0;
    };
    VirtualPolyhedronProjections vp(E, nu, vf, dofs, mesh, funcx, funcy, funcz, k);

    SolverVEM solvervem(k, 3 * DOFS.getnumDofs(), dofs, vp);

    // std::cout<<solvervem.getRightHandSide()<<std::endl;

    std::vector<std::function<real(real, real, real)>> HomogeneousDirichletBC;
    HomogeneousDirichletBC.emplace_back([](real x, real y, real z)
                                        { return z; });
    /*
    HomogeneousDirichletBC.emplace_back([](real x, real y, real z)
                                        { return x; });
    HomogeneousDirichletBC.emplace_back([](real x, real y, real z)
                                        { return y; });
    HomogeneousDirichletBC.emplace_back([](real x, real y, real z)
                                        { return z; });
    HomogeneousDirichletBC.emplace_back([](real x, real y, real z)
                                        { return x - 1.0; });
    HomogeneousDirichletBC.emplace_back([](real x, real y, real z)
                                        { return y - 1.0; });
    HomogeneousDirichletBC.emplace_back([](real x, real y, real z)
                                        { return z - 1.0; });
    */
    solvervem.enforceHomogeneousDirichletBC(mesh, DOFS, HomogeneousDirichletBC);

    /*
        Eigen::MatrixXi A = Eigen::MatrixXi::Random(4, 6);
        std::cout << "Initial matrix A:\n"
                  << A << "\n\n";
        std::cout << "A(all,{4,2,5,5,3}):\n"
                  << A(all, {4, 2, 5, 5, 3}) << "\n\n";
    */

    solvervem.solve();
    //auto U=solvervem.getSolutionDisplacementsEig();
    std::vector<real> sol=solvervem.getSolutionDisplacements();
    
    std::cout<<"PRINTING UVECTOR"<<std::endl;
    std::cout<<sol.size()<<std::endl;
    std::cout<<sol[0]<<std::endl;
    for (const auto& s:sol)
    {
        std::cout<<s<<std::endl;
    }
    

    return 0;
}