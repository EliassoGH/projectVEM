#include "problem.hpp"

Problem::Problem(const char *parametersFilename, bool printError)
{

    // READ PARAMETERS
    std::cout << "Reading parameters..." << std::endl;
    Parameters parameters(parametersFilename);
    parameters.print();
    std::cout << std::endl;
    const unsigned int order = parameters.getOrder();

    // READ MESH
    std::string meshName = parameters.getInputMesh();
    std::cout << "Reading mesh " << meshName << "..." << std::endl;
    Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> mesh(meshName);
    mesh.print();
    std::cout << std::endl;

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

    // SOLVE THE SYSTEM
    std::cout << "Solving the system..." << std::endl;
    solvervem.solve();

    // EXPORT SOLUTION
    std::cout << "Writing output file..." << std::endl;
    plotting::export_results(mesh, dofs, PProj, solvervem.getSolutionDisplacements());
    std::cout<<std::endl;

    if (printError)
    {
        std::cout << std::setprecision(16);
        std::cout << "L2 strain error = " << solvervem.computeStrainError(mesh, dofs, PProj, parameters.getExactStrains()) << std::endl;
    }
}