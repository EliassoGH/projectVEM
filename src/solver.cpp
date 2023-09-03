#include "solver.hpp"

real computeSparseMatrixTrace(const Eigen::SparseMatrix<real> &matrix)
{
    real trace = 0.0;
    for (Eigen::Index i = 0; i < matrix.rows(); ++i)
    {
        trace += matrix.coeff(i, i);
    }
    return trace;
}

void Solver::solve()
{
    K.prune(1e-12, 1);

    K.finalize();
    K.makeCompressed();

    Eigen::ConjugateGradient<Eigen::SparseMatrix<real>, Eigen::Upper, Eigen::IncompleteCholesky<real>> solver;
    solver.compute(K);

    if (solver.info() != Eigen::Success)
    {
        throw std::runtime_error("Solver factorize failed");
    }

    U = solver.solve(F);

    if (solver.info() != Eigen::Success)
    {
        throw std::runtime_error("Solver solving failed");
    }

}

const Eigen::VectorXd &Solver::getRightHandSide() const
{
    return F;
}

const Eigen::VectorXd &Solver::getSolutionDisplacementsEig() const
{
    return U;
}

std::vector<real> Solver::getSolutionDisplacements() const
{
    std::vector<real> U_vec;
    U_vec.resize(U.size());
    Eigen::Map<Eigen::VectorXd>(U_vec.data(), U_vec.size()) = U;
    return U_vec;
}

real Solver::computeH1error(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
                            const std::function<std::array<real, 3>(real, real, real)> &Uex_func)
{
    Eigen::VectorXd U_ex(3 * mesh.numVertices());
    for (std::size_t i = 0; i < mesh.numVertices(); i++)
    {
        auto v = mesh.getVertex(i);
        auto Ui = Uex_func(v[0], v[1], v[2]);
        U_ex[3 * i] = Ui[0];
        U_ex[3 * i + 1] = Ui[0];
        U_ex[3 * i + 2] = Ui[0];
    }

    return (U_ex - U.head(3 * mesh.numVertices())).transpose() * K.block(0, 0, 3 * mesh.numVertices(), 3 * mesh.numVertices()) * (U_ex - U.head(3 * mesh.numVertices()));
}

SolverVEM::SolverVEM(const Parameters &parameters,
                     const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
                     const VirtualDofsCollection &DOFS,
                     const LocalVirtualDofsCollection &dofs,
                     const VirtualPolyhedronProjections &vp)
{
    enforceHomogeneousDirichletBC(mesh, DOFS, parameters.getHomogeneousDirichletBC());
    std::size_t system_dim = 3 * DOFS.getnumDofs();
    unsigned int order = parameters.getOrder();
    K.resize(system_dim, system_dim);
    F.setZero(system_dim);

    std::vector<Eigen::Triplet<real>> tripletList_glob;
#pragma omp declare reduction(merge : std::vector<Eigen::Triplet<real>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
    std::vector<Eigen::Triplet<real>> tripletList;
// #pragma omp parallel
//{
// std::vector<Eigen::Triplet<real>> tripletList;
//#pragma omp parallel for
#pragma omp parallel for reduction(merge: tripletList)
    for (std::size_t p = 0; p < dofs.numLocalDofsCollection(); p++)
    {
        const Eigen::SparseMatrix<real> &C = vp.getPolyhedronProjections().at(p);
        const Eigen::SparseMatrix<real> &E = vp.getElasticMatrices().at(p);
        const Eigen::SparseMatrix<real> &Tdr = vp.getDeformationRBMMatrices().at(p);
        const Eigen::VectorXd &f = vp.getforcingProjections().at(p);
        // Assemble K_ec
        Eigen::SparseMatrix<real> K_ec = C.transpose() * E.selfadjointView<Eigen::Upper>() * C;

        // Transform to dense (Tdr^T)*Tdr so it can be inverted
        Eigen::MatrixXd TdrTTdr = Tdr.transpose() * Tdr;
        // Create identity matrix
        Eigen::SparseMatrix<real> I(K_ec.rows(), K_ec.rows());
        I.setIdentity();

        // Assemble K_es
        Eigen::MatrixXd K_es = 0.5 * computeSparseMatrixTrace(K_ec) * (I - (Tdr * TdrTTdr.inverse() * Tdr.transpose()));
        // Assemble K_e
        Eigen::MatrixXd K_e = K_es;
        K_e += K_ec;

        auto pdofs = dofs.getLocalDofs(p);
        const auto npdofs = pdofs.getnumDofs();
        for (std::size_t i = 0; i < npdofs; i++)
        {
            std::size_t IDi = pdofs.getID(i) * 3;
            if (isConstrained[IDi])
            {
                for (std::size_t ii = 0; ii < 3; ii++)
                {
                    // diagonal of diagonal blocks
                    tripletList.push_back(Eigen::Triplet<real>(IDi + ii, IDi + ii, 1.0));
                }
                continue;
            }
            std::size_t idi = 3 * i;

            for (std::size_t ii = 0; ii < 3; ii++)
            {
                // diagonal of diagonal blocks
                tripletList.push_back(Eigen::Triplet<real>(IDi + ii, IDi + ii, K_e(idi + ii, idi + ii)));
            }
            // extradiagonals of diagonal blocks
            tripletList.push_back(Eigen::Triplet<real>(IDi, IDi + 1, K_e(idi, idi + 1)));
            tripletList.push_back(Eigen::Triplet<real>(IDi, IDi + 2, K_e(idi, idi + 2)));
            tripletList.push_back(Eigen::Triplet<real>(IDi + 1, IDi + 2, K_e(idi + 1, idi + 2)));

            for (std::size_t j = i + 1; j < npdofs; j++)
            {
                std::size_t IDj = pdofs.getID(j) * 3;
                if (isConstrained[IDj])
                {
                    continue;
                }
                std::size_t idj = 3 * j;
                if (IDj > IDi)
                {
                    for (std::size_t ii = 0; ii < 3; ii++)
                    {
                        // diagonal of extradiagonal blocks
                        tripletList.push_back(Eigen::Triplet<real>(IDi + ii, IDj + ii, K_e(idi + ii, idj + ii)));
                    }
                    // extradiagonals of extradiagonal blocks
                    tripletList.push_back(Eigen::Triplet<real>(IDi, IDj + 1, K_e(idi, idj + 1)));
                    tripletList.push_back(Eigen::Triplet<real>(IDi, IDj + 2, K_e(idi, idj + 2)));
                    tripletList.push_back(Eigen::Triplet<real>(IDi + 1, IDj + 2, K_e(idi + 1, idj + 2)));
                    tripletList.push_back(Eigen::Triplet<real>(IDi + 1, IDj, K_e(idi + 1, idj)));
                    tripletList.push_back(Eigen::Triplet<real>(IDi + 2, IDj, K_e(idi + 2, idj)));
                    tripletList.push_back(Eigen::Triplet<real>(IDi + 2, IDj + 1, K_e(idi + 2, idj + 1)));
                }
                else
                {
                    for (std::size_t ii = 0; ii < 3; ii++)
                    {
                        // diagonal of extradiagonal blocks
                        tripletList.push_back(Eigen::Triplet<real>(IDj + ii, IDi + ii, K_e(idj + ii, idi + ii)));
                    }
                    // extradiagonals of extradiagonal blocks
                    tripletList.push_back(Eigen::Triplet<real>(IDj, IDi + 1, K_e(idj, idi + 1)));
                    tripletList.push_back(Eigen::Triplet<real>(IDj, IDi + 2, K_e(idj, idi + 2)));
                    tripletList.push_back(Eigen::Triplet<real>(IDj + 1, IDi + 2, K_e(idj + 1, idi + 2)));
                    tripletList.push_back(Eigen::Triplet<real>(IDj + 1, IDi, K_e(idj + 1, idi)));
                    tripletList.push_back(Eigen::Triplet<real>(IDj + 2, IDi, K_e(idj + 2, idi)));
                    tripletList.push_back(Eigen::Triplet<real>(IDj + 2, IDi + 1, K_e(idj + 2, idi + 1)));
                }
            }
#pragma omp critical
            {
                if (order == 1)
                {
                    F[IDi] += f[0];
                    F[IDi + 1] += f[1];
                    F[IDi + 2] += f[2];
                }
            }
        }
#pragma omp critical
        {
            for (std::size_t i = 0; i < pdofs.getnumPdofs(); i++)
            {
                F[3 * pdofs.getID(npdofs - pdofs.getnumPdofs() + i)] += f[i];
                F[3 * pdofs.getID(npdofs - pdofs.getnumPdofs() + i) + 1] += f[i + 1];
                F[3 * pdofs.getID(npdofs - pdofs.getnumPdofs() + i) + 2] += f[i + 2];
            }
        }
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());
}

void SolverVEM::enforceHomogeneousDirichletBC(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
                                              const VirtualDofsCollection &DOFS,
                                              const std::vector<std::function<real(real, real, real)>> &constraintFunctions)
{
    const auto nVDOFS = DOFS.getnumVdofs();
    const auto nEDOFS = DOFS.getnumEdofs();
    const auto nFDOFS = DOFS.getnumFdofs();
    //const auto nPDOFS = DOFS.getnumPdofs();
    isConstrained.resize(3 * DOFS.getnumDofs(), false);

    bool unconstrained = true;
    for (std::size_t j = 0; j < nVDOFS; ++j)
    {
        unconstrained = true;
        // Iterate through each constraint function
        for (std::size_t i = 0; i < constraintFunctions.size(); ++i)
        {
            const auto &constraintFunc = constraintFunctions[i];
            const auto &VDOF = DOFS.getDof<VertexDof>(j)->getVertex();
            real constraintValue = constraintFunc(VDOF[0], VDOF[1], VDOF[2]);

            if (std::abs(constraintValue) < tolerance)
            {
                unconstrained = false;
                // add the dof id to the constrained dofs vector
                /*
                constrainedDofs.push_back(3 * j);
                constrainedDofs.push_back(3 * j + 1);
                constrainedDofs.push_back(3 * j + 2);
                */
                // add value 0.0 to the constrained dof values vector
                /*
                constrainedDofsValues.push_back(0.0);
                constrainedDofsValues.push_back(0.0);
                constrainedDofsValues.push_back(0.0);
                */
                // add true to logical isConstrained
                isConstrained[3 * j] = true;
                isConstrained[3 * j + 1] = true;
                isConstrained[3 * j + 2] = true;
                break;
            }
        }
        if (unconstrained)
        {
            // add the dof id to the unconstrained dofs vector
            /*
            unconstrainedDofs.push_back(3 * j);
            unconstrainedDofs.push_back(3 * j + 1);
            unconstrainedDofs.push_back(3 * j + 2);
            */
        }
    }
    for (std::size_t e = 0; e < nEDOFS; ++e)
    {
        unconstrained = true;
        auto j = nVDOFS + e;
        const auto &EDOF = DOFS.getDof<EdgeDof>(j)->getGaussLobattoPoint();
        // Iterate through each constraint function
        for (std::size_t i = 0; i < constraintFunctions.size(); ++i)
        {
            const auto &constraintFunc = constraintFunctions[i];
            real constraintValue = constraintFunc(EDOF[0], EDOF[1], EDOF[2]);

            if (std::abs(constraintValue) < tolerance)
            {
                unconstrained = false;
                // add the dof id to the constrained dofs vector
                /*
                constrainedDofs.push_back(3 * j);
                constrainedDofs.push_back(3 * j + 1);
                constrainedDofs.push_back(3 * j + 2);
                */
                // add value 0.0 to the constrained dof values vector
                /*
                constrainedDofsValues.push_back(0.0);
                constrainedDofsValues.push_back(0.0);
                constrainedDofsValues.push_back(0.0);
                */
                // add true to logical isConstrained
                isConstrained[3 * j] = true;
                isConstrained[3 * j + 1] = true;
                isConstrained[3 * j + 2] = true;
                break;
            }
        }
        if (unconstrained)
        {
            // add the dof id to the unconstrained dofs vector
            /*
            unconstrainedDofs.push_back(3 * j);
            unconstrainedDofs.push_back(3 * j + 1);
            unconstrainedDofs.push_back(3 * j + 2);
            */
        }
    }
    for (std::size_t f = 0; f < nFDOFS; ++f)
    {
        unconstrained = true;
        auto j = nVDOFS + nEDOFS + f;
        auto FID = DOFS.getDof<FaceDof>(j)->getId();
        const auto &v1 = mesh.getPolygon(FID)[0][0];
        const auto &v2 = mesh.getPolygon(FID)[0][1];
        const auto &v3 = mesh.getPolygon(FID)[1][1];

        // Iterate through each constraint function
        for (std::size_t i = 0; i < constraintFunctions.size(); ++i)
        {
            const auto &constraintFunc = constraintFunctions[i];
            real constraintValuev1 = constraintFunc(v1[0], v1[1], v1[2]);
            real constraintValuev2 = constraintFunc(v2[0], v2[1], v2[2]);
            real constraintValuev3 = constraintFunc(v3[0], v3[1], v3[2]);

            if ((std::abs(constraintValuev1) < tolerance) &&
                (std::abs(constraintValuev2) < tolerance) &&
                (std::abs(constraintValuev3) < tolerance))
            {
                unconstrained = false;
                // add the dof id to the constrained dofs vector
                /*
                constrainedDofs.push_back(3 * j);
                constrainedDofs.push_back(3 * j + 1);
                constrainedDofs.push_back(3 * j + 2);
                */
                // add value 0.0 to the constrained dof values vector
                /*
                constrainedDofsValues.push_back(0.0);
                constrainedDofsValues.push_back(0.0);
                constrainedDofsValues.push_back(0.0);
                */
                // add true to logical isConstrained
                isConstrained[3 * j] = true;
                isConstrained[3 * j + 1] = true;
                isConstrained[3 * j + 2] = true;
                break;
            }
        }
        if (unconstrained)
        {
            // add the dof id to the unconstrained dofs vector
            /*
            unconstrainedDofs.push_back(3 * j);
            unconstrainedDofs.push_back(3 * j + 1);
            unconstrainedDofs.push_back(3 * j + 2);
            */
        }
    }
    /*
    for (std::size_t p = 0; p < nPDOFS; p++)
    {
        auto j = nVDOFS + nEDOFS + nFDOFS + p;
        // add the dof id to the unconstrained dofs vector
        
        unconstrainedDofs.push_back(3 * j);
        unconstrainedDofs.push_back(3 * j + 1);
        unconstrainedDofs.push_back(3 * j + 2);
        
    }
    */
}

real SolverVEM::computeStrainError(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
                                   const LocalVirtualDofsCollection &dofs,
                                   const VirtualPolyhedronProjections &vp,
                                   const std::function<std::array<real, 6>(real, real, real)> &EpsEx_func)
{
    real error = 0.0;
    unsigned int order = vp.getOrder();
    const auto m_kminus1 = Monomial3D::getMonomialsOrdered(order - 1);
    // For every polyhedron
    for (std::size_t p = 0; p < dofs.numLocalDofsCollection(); p++)
    {
        auto P = mesh.getPolyhedron(p);
        Point3D X_P = getPolyhedronCentroid(P);
        real h_P = P.getDiameter();
        const Eigen::SparseMatrix<real> &C = vp.getPolyhedronProjections().at(p);
        auto pdofs = dofs.getLocalDofs(p);
        Eigen::VectorXd U_loc(3 * pdofs.getnumDofs());
        // Find strains coefficients
        for (std::size_t i = 0; i < pdofs.getnumDofs(); i++)
        {
            U_loc[3 * i] = U[3 * pdofs.getID(i)];
            U_loc[3 * i + 1] = U[3 * pdofs.getID(i) + 1];
            U_loc[3 * i + 2] = U[3 * pdofs.getID(i) + 2];
        }
        Eigen::VectorXd EpsApprox_coeff = C * U_loc;

        auto epsErr_func = [&EpsEx_func, &EpsApprox_coeff, &m_kminus1, &X_P, &h_P](real x, real y, real z)
        {
            real result = 0.0;
            // For each of the 6 components of the strain field
            for (std::size_t j = 0; j < 6; j++)
            {
                real EpsApprox = 0.0;
                // For every monomial of the strain field
                for (std::size_t i = 0; i < m_kminus1.size(); i++)
                {
                    EpsApprox += EpsApprox_coeff[6 * i + j] * m_kminus1[i].evaluate((Point3D(x, y, z) - X_P) / h_P);
                }
                result += std::pow((EpsEx_func(x, y, z)[j] - EpsApprox), 2);
            };
            return result;
        };
        error += (integrateFunctionOverPolyhedron(epsErr_func, P, order + 1));
    }
    return std::sqrt(error);
}