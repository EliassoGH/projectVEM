#ifndef __PROBLEM_HPP_
#define __PROBLEM_HPP_

#include "mesh.hpp"
#include "parameters.hpp"
#include "monomial.hpp"
#include "virtualDofs.hpp"
#include "virtualProjections.hpp"
#include "solver.hpp"
#include "export_results.hpp"
#include <iostream>
#include <chrono>

class Problem
{
public:
    Problem(const char *meshFilename = "mesh.geo", const char *parametersFilename = "parameters.dat")
    {
        // READ MESH
        std::cout << "Reading mesh..." << std::endl;
        Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> mesh(meshFilename);
        mesh.print();
        std::cout << std::endl;

        // READ PARAMETERS
        std::cout << "Reading parameters..." << std::endl;
        Parameters parameters(parametersFilename);
        parameters.print();
        std::cout << std::endl;
        const unsigned int order = parameters.getOrder();
        if (parameters.getInitializeFaceIntegrals())
            MonomialsFaceIntegralsCache::initialize(mesh, 2 * order - 1);

        // CREATE DEGREES OF FREEDOM
        std::cout << "Creating global degrees of freedom..." << std::endl;
        VirtualDofsCollection DOFS(mesh, order);
        DOFS.print();
        std::cout << std::endl;
        std::cout << "Assigning local degrees of freedom..." << std::endl;
        LocalVirtualDofsCollection dofs(mesh, DOFS);

        // COMPUTE VIRTUAL PROJECTIONS
        std::cout << "Computing face projections..." << std::endl;
        VirtualFaceProjections FProj(DOFS, mesh, order);
        std::cout << "Computing polyhedron projections..." << std::endl;
        VirtualPolyhedronProjections PProj(parameters, FProj, dofs, mesh);

        // ASSEMBLE THE SYSTEM
        std::cout << "Assembling the system..." << std::endl;
        SolverVEM solvervem(parameters, mesh, DOFS, dofs, PProj);
/*
        // ENFORCE HOMOGENEOUS DIRICHLET BOUNDARY CONDITIONS
        std::cout << "Enforcing homogeneous dirichlet boundary conditions..." << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        solvervem.enforceHomogeneousDirichletBC(mesh, DOFS, parameters.getHomogeneousDirichletBC());
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "time taken: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
*/
        // SOLVE THE SYSTEM
        std::cout << "Solving the system..." << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        solvervem.solve();
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "time taken: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;

        // GET SOLUTION
        std::cout << "Retrieving solution..." << std::endl;
        start = std::chrono::high_resolution_clock::now();
        std::vector<real> solution = solvervem.getSolutionDisplacements();
        // std::cout<<sol<<std::endl;
        end = std::chrono::high_resolution_clock::now();
        std::cout << "time taken: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;

        // EXPORT SOLUTION
        std::cout << "Writing output file..." << std::endl;
        start = std::chrono::high_resolution_clock::now();
        plotting::export_results(mesh, solution);
        end = std::chrono::high_resolution_clock::now();
        std::cout << "time taken: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;

        start = std::chrono::high_resolution_clock::now();
        std::cout << "error" << solvervem.computeH1error(mesh, parameters.getExactSolution())<<std::endl;
        end = std::chrono::high_resolution_clock::now();
        std::cout << "time taken: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
    }
};

class ConvergenceAnalysis
{
protected:
    std::vector<real> mesh_sizes;
    std::vector<real> method_orders;
    std::vector<real> errors;

public:
    ConvergenceAnalysis() {}
};

class Convergence_h_refinement : public ConvergenceAnalysis
{
public:
    Convergence_h_refinement()
    {
        // READ PARAMETERS
        std::cout << "Reading parameters..." << std::endl;
        Parameters parameters("parameters.dat");
        parameters.print();
        std::cout << std::endl;
        const unsigned int order = parameters.getOrder();

        std::string meshName = "N2.geo";
        // READ MESH
        std::cout << "Reading mesh " << meshName << "..." << std::endl;
        Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> mesh(meshName);
        mesh.print();
        std::cout << std::endl;

        if (parameters.getInitializeFaceIntegrals())
            MonomialsFaceIntegralsCache::initialize(mesh, 2 * order - 1);

        std::cout << "Solving..." << std::endl;
        // CREATE DEGREES OF FREEDOM
        VirtualDofsCollection DOFS(mesh, order);
        LocalVirtualDofsCollection dofs(mesh, DOFS);

        // COMPUTE VIRTUAL PROJECTIONS
        VirtualFaceProjections FProj(DOFS, mesh, order);
        VirtualPolyhedronProjections PProj(parameters, FProj, dofs, mesh);

        // ASSEMBLE THE SYSTEM
        //SolverVEM solvervem(order, 3 * DOFS.getnumDofs(), dofs, PProj);
        SolverVEM solvervem(parameters, mesh, DOFS, dofs, PProj);

        // ENFORCE HOMOGENEOUS DIRICHLET BOUNDARY CONDITIONS
        //solvervem.enforceHomogeneousDirichletBC(mesh, DOFS, parameters.getHomogeneousDirichletBC());

        // SOLVE THE SYSTEM
        solvervem.solve();

        // GET SOLUTION
        //std::vector<real> solution = solvervem.getSolutionDisplacements();
        //std::cout<<solvervem.getSolutionDisplacementsEig()<<std::endl;

        // COMPUTE ERROR
        /*
        errors.push_back(solvervem.computeH1error(mesh, parameters.getExactSolution()));
        mesh_sizes.push_back(mesh.getSize());
        for (std::size_t i = 0; i < mesh_sizes.size(); i++)
        {
            std::cout << std::setprecision(10);
            std::cout << mesh_sizes[i] << " " << errors[i] << std::endl;
        }
        */
        std::cout<<std::setprecision(16);
        std::cout<<"error is "<<solvervem.computeStrainError(mesh,dofs,PProj,parameters.getExactStrains())<<std::endl;
    }
};

#endif // __PROBLEM_HPP_