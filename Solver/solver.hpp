#ifndef __SOLVER_HPP_
#define __SOLVER_HPP_

#include "mesh.hpp"
#include "monomial.hpp"
#include "integration.hpp"
#include "virtualDofs.hpp"
#include "virtualProjections.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <omp.h>
#include <map>
#include <vector>
#include <functional>
#include <chrono>
#include <Eigen/IterativeLinearSolvers>

using namespace IntegrationMonomial;
using namespace Gauss;

real computeSparseMatrixTrace(const Eigen::SparseMatrix<real> &matrix)
{
    real trace = 0.0;
    for (Eigen::Index i = 0; i < matrix.rows(); ++i)
    {
        trace += matrix.coeff(i, i);
    }
    return trace;
}

template <typename T, typename T2>
Eigen::Matrix<typename T2::Scalar, T::RowsAtCompileTime, T::ColsAtCompileTime, T::Options>
extract(const Eigen::DenseBase<T2> &full, const Eigen::DenseBase<T> &ind)
{
    using target_t = Eigen::Matrix<typename T2::Scalar, T::RowsAtCompileTime, T::ColsAtCompileTime, T::Options>;
    int num_indices = ind.innerSize();
    target_t target(num_indices);
    for (int i = 0; i < num_indices; i++)
    {
        target[i] = full[ind[i]];
    }
    return target;
}

class Solver
{
protected:
    Eigen::SparseMatrix<real> K;                // System matrix, upper part
    Eigen::VectorXd F;                          // Right-hand side vector
    std::vector<std::size_t> unconstrainedDofs; // vector storing uncosntrained dofs
    std::vector<std::size_t> constrainedDofs;   // vector storing constrained dofs
    std::vector<real> constrainedDofsValues;    // vector storing constrained dofs
    Eigen::VectorXd U;                          // solution
    std::vector<bool> isConstrained;            // vector storing true if dof is constrained

public:
    Solver() {}
    Solver(const Eigen::SparseMatrix<real> &K, const Eigen::VectorXd &F) : K(K), F(F) {}

    // Eigen::VectorXd solve() const
    void solve()
    {
        /*
        Eigen::Map<const Eigen::VectorXd> constrainedDofsValuesEig(constrainedDofsValues.data(), constrainedDofsValues.size());
        Eigen::MatrixXd K_dense = K.toDense();
        K_dense = (m.array() >= 5).select(K_dense, m);
        Eigen::MatrixXd K_reduced=K_dense(Eigen::all,constrainedDofsValuesEig);
        std::cout<<K_reduced<<std::endl;
        */

        // Eigen::MatrixXd K_dense = K.toDense();
        //  std::cout<<K<<std::endl<<std::endl;
        // std::cout << K_dense << std::endl<< std::endl;
        //  std::cout<<F<<std::endl<<std::endl;
        //  std::cout<<"constrained dofs size "<<constrainedDofs.size()<<std::endl;

        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();

        unsigned int elapsed_time = 0;
        /*
        for (std::size_t i = 0; i < constrainedDofs.size(); i++)
        {
            std::size_t constrainedDof = constrainedDofs[i];
            K.prune([&constrainedDof](std::size_t i, std::size_t j, real)
                    { return i != constrainedDof && j != constrainedDof; });
            auto start = std::chrono::high_resolution_clock::now();
            K.insert(constrainedDof, constrainedDof) = 1.0;
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            // std::cout<<"constrained value "<<constrainedDofsValues[i]<<" at "<<constrainedDof;
            F[constrainedDof] = constrainedDofsValues[i];
            // std::cout<<", f inserted "<<F[constrainedDof]<<std::endl;
        }
        */
        std::cout << "modifying matrix for BC: " << elapsed_time << std::endl;

        // Eigen::MatrixXd K_dense1 = K.toDense();
        //  std::cout<<K_dense1<<std::endl<<std::endl;
        //  std::cout<<F<<std::endl<<std::endl;

        /*
        for (std::size_t i=0;i<constrainedDofs.size();i++)
        {
            std::cout<<constrainedDofs[i]<<" "<<constrainedDofsValues[i]<<std::endl;
        }
        */

        // auto F_reduced=F-K_dense(constrainedDofs,constrainedDofs)*constrainedDofsValuesEig;
        // std::cout<<F_reduced<<std::endl;

        start = std::chrono::high_resolution_clock::now();
        K.pruned();
        K.finalize();
        K.makeCompressed();
        // Eigen::SparseLU<Eigen::SparseMatrix<real>> solver;
        Eigen::ConjugateGradient<Eigen::SparseMatrix<real>, Eigen::Upper> solver;
        solver.compute(K);
        // std::cout << "factorized" << std::endl;

        if (solver.info() != Eigen::Success)
        {
            throw std::runtime_error("Solver factorize failed");
        }
        // std::cout << "about to solve" << std::endl;
        U = solver.solve(F);

        if (solver.info() != Eigen::Success)
        {
            throw std::runtime_error("Solver solving failed");
        }
        end = std::chrono::high_resolution_clock::now();
        std::cout << "actual solving: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
    }

    const Eigen::VectorXd &getRightHandSide() const
    {
        return F;
    }

    const Eigen::VectorXd &getSolutionDisplacementsEig() const
    {
        return U;
    }

    std::vector<real> getSolutionDisplacements() const
    {
        std::vector<real> U_vec;
        U_vec.resize(U.size());
        Eigen::Map<Eigen::VectorXd>(U_vec.data(), U_vec.size()) = U;
        return U_vec;
        // return std::vector<real>(&U[0], U.data() + U.size());
        // return (std::vector<real>(U.data(), U.data() + U.size()));
    }

    real computeH1error(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
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

        // Eigen::SparseMatrix<real> Kv=K.block(0,0,3*mesh.numVertices(),3*mesh.numVertices());
        /*
        std::cout<<"approx solution at vertices"<<std::endl;
        std::cout<<U.head(3*mesh.numVertices())<<std::endl;
        std::cout<<"exact  solution at vertices"<<std::endl;
        std::cout<<U_ex<<std::endl;
        std::cout<<"difference"<<std::endl;
        std::cout<<(U_ex-U.head(3*mesh.numVertices()))<<std::endl;
        */
        return (U_ex - U.head(3 * mesh.numVertices())).transpose() * K.block(0, 0, 3 * mesh.numVertices(), 3 * mesh.numVertices()) * (U_ex - U.head(3 * mesh.numVertices()));
    }
};

class SolverVEM : public Solver
{
private:
    static constexpr real tolerance = 1e-12;

public:
    SolverVEM(const Parameters &parameters,
              const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
              const VirtualDofsCollection &DOFS,
              const LocalVirtualDofsCollection &dofs,
              const VirtualPolyhedronProjections &vp)
    {
        enforceHomogeneousDirichletBC(mesh, DOFS, parameters.getHomogeneousDirichletBC());
        std::size_t system_dim = 3 * DOFS.getnumDofs();
        unsigned int order = parameters.getOrder();
        std::vector<Eigen::Triplet<real>> tripletList;
        // tripletList.reserve(system_dim * system_dim / 2);
        tripletList.reserve(system_dim);
        K.resize(system_dim, system_dim);
        F.setZero(system_dim);
        unsigned int compute_K_e = 0.0;
        unsigned int place_K_e = 0.0;

        // std::cout << "HERE1" << std::endl;
        for (std::size_t p = 0; p < dofs.numLocalDofsCollection(); p++)
        {
            auto start = std::chrono::high_resolution_clock::now();
            /*
            std::map<std::size_t, Eigen::SparseMatrix<real>> polyhedronProjections;    // stores matrices C
            std::map<std::size_t, Eigen::SparseMatrix<real>> elastic_matrices;         // stores matrices E
            std::map<std::size_t, Eigen::SparseMatrix<real>> deformation_RBM_matrices; // stores matrices Tdr
            std::map<std::size_t, Eigen::VectorXd> forcingProjections;                 // stores polynomials b_hat
            */
            const Eigen::SparseMatrix<real> &C = vp.getPolyhedronProjections().at(p);
            const Eigen::SparseMatrix<real> &E = vp.getElasticMatrices().at(p);
            const Eigen::SparseMatrix<real> &Tdr = vp.getDeformationRBMMatrices().at(p);
            const Eigen::VectorXd &f = vp.getforcingProjections().at(p);
            Eigen::SparseMatrix<real> K_ec = C.transpose() * E.selfadjointView<Eigen::Upper>() * C;

            // Create identity matrix
            auto I = Eigen::MatrixXd::Identity(K_ec.rows(), K_ec.cols());
            // Transform to dense (Tdr^T)*Tdr so it can be inverted
            Eigen::MatrixXd TdrTTrd = Tdr.transpose() * Tdr;
            // Assemble K_es
            Eigen::MatrixXd K_es = 0.5 * computeSparseMatrixTrace(K_ec) * (I - (Tdr * TdrTTrd.inverse() * Tdr.transpose()));
            Eigen::MatrixXd K_e = K_es;
            K_e += K_ec;
            // Eigen::SparseMatrix<real> K_e_sparse = K_e.sparseView();
            //  std::cout<<K_ec.toDense()<<std::endl<<std::endl;
            //  std::cout<<K_es<<std::endl<<std::endl;
            //  std::cout << K_e << std::endl;

            auto end = std::chrono::high_resolution_clock::now();
            // std::cout << "compute K_e: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
            compute_K_e += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            start = std::chrono::high_resolution_clock::now();

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
                /*
                tripletList.push_back(Eigen::Triplet<real>(IDi + 1, IDi, K_e(idi + 1, idi)));
                tripletList.push_back(Eigen::Triplet<real>(IDi + 2, IDi, K_e(idi + 2, idi)));
                tripletList.push_back(Eigen::Triplet<real>(IDi + 2, IDi + 1, K_e(idi + 2, idi + 1)));
                */

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
                            // tripletList.push_back(Eigen::Triplet<real>(IDj + ii, IDi + ii, K_e(idj + ii, idi + ii)));
                        }
                        // extradiagonals of extradiagonal blocks
                        tripletList.push_back(Eigen::Triplet<real>(IDi, IDj + 1, K_e(idi, idj + 1)));
                        tripletList.push_back(Eigen::Triplet<real>(IDi, IDj + 2, K_e(idi, idj + 2)));
                        tripletList.push_back(Eigen::Triplet<real>(IDi + 1, IDj + 2, K_e(idi + 1, idj + 2)));
                        tripletList.push_back(Eigen::Triplet<real>(IDi + 1, IDj, K_e(idi + 1, idj)));
                        tripletList.push_back(Eigen::Triplet<real>(IDi + 2, IDj, K_e(idi + 2, idj)));
                        tripletList.push_back(Eigen::Triplet<real>(IDi + 2, IDj + 1, K_e(idi + 2, idj + 1)));
                        /*
                        tripletList.push_back(Eigen::Triplet<real>(IDj, IDi + 1, K_e(idj, idi + 1)));
                        tripletList.push_back(Eigen::Triplet<real>(IDj, IDi + 2, K_e(idj, idi + 2)));
                        tripletList.push_back(Eigen::Triplet<real>(IDj + 1, IDi + 2, K_e(idj + 1, idi + 2)));
                        tripletList.push_back(Eigen::Triplet<real>(IDj + 1, IDi, K_e(idj + 1, idi)));
                        tripletList.push_back(Eigen::Triplet<real>(IDj + 2, IDi, K_e(idj + 2, idi)));
                        tripletList.push_back(Eigen::Triplet<real>(IDj + 2, IDi + 1, K_e(idj + 2, idi + 1)));
                        */
                    }
                    else
                    {
                        for (std::size_t ii = 0; ii < 3; ii++)
                        {
                            // diagonal of extradiagonal blocks
                            // tripletList.push_back(Eigen::Triplet<real>(IDi + ii, IDj + ii, K_e(idi + ii, idj + ii)));
                            tripletList.push_back(Eigen::Triplet<real>(IDj + ii, IDi + ii, K_e(idj + ii, idi + ii)));
                        }
                        // extradiagonals of extradiagonal blocks
                        /*
                        tripletList.push_back(Eigen::Triplet<real>(IDi, IDj + 1, K_e(idi, idj + 1)));
                        tripletList.push_back(Eigen::Triplet<real>(IDi, IDj + 2, K_e(idi, idj + 2)));
                        tripletList.push_back(Eigen::Triplet<real>(IDi + 1, IDj + 2, K_e(idi + 1, idj + 2)));
                        tripletList.push_back(Eigen::Triplet<real>(IDi + 1, IDj, K_e(idi + 1, idj)));
                        tripletList.push_back(Eigen::Triplet<real>(IDi + 2, IDj, K_e(idi + 2, idj)));
                        tripletList.push_back(Eigen::Triplet<real>(IDi + 2, IDj + 1, K_e(idi + 2, idj + 1)));
                        */
                        tripletList.push_back(Eigen::Triplet<real>(IDj, IDi + 1, K_e(idj, idi + 1)));
                        tripletList.push_back(Eigen::Triplet<real>(IDj, IDi + 2, K_e(idj, idi + 2)));
                        tripletList.push_back(Eigen::Triplet<real>(IDj + 1, IDi + 2, K_e(idj + 1, idi + 2)));
                        tripletList.push_back(Eigen::Triplet<real>(IDj + 1, IDi, K_e(idj + 1, idi)));
                        tripletList.push_back(Eigen::Triplet<real>(IDj + 2, IDi, K_e(idj + 2, idi)));
                        tripletList.push_back(Eigen::Triplet<real>(IDj + 2, IDi + 1, K_e(idj + 2, idi + 1)));
                    }
                }

                if (order == 1)
                {
                    F[IDi] += f[0];
                    F[IDi + 1] += f[1];
                    F[IDi + 2] += f[2];
                }
            }

            for (std::size_t i = 0; i < pdofs.getnumPdofs(); i++)
            {
                F[3 * pdofs.getID(npdofs - pdofs.getnumPdofs() + i)] += f[i];
                F[3 * pdofs.getID(npdofs - pdofs.getnumPdofs() + i) + 1] += f[i + 1];
                F[3 * pdofs.getID(npdofs - pdofs.getnumPdofs() + i) + 2] += f[i + 2];
            }
            K.setFromTriplets(tripletList.begin(), tripletList.end());
            /*
            for (std::size_t i = 0; i < npdofs; i++)
            {
                std::size_t IDi = pdofs.getID(i) * 3;
                std::size_t idi = 3 * i;

                for (std::size_t ii = 0; ii < 3; ii++)
                {
                    // std::cout << "inserting local (" << 3 * i + ii << ", " << 3 * i + ii << ") into global (" << 3 * pdofs.getID(i) + ii << ", " << 3 * pdofs.getID(i) + ii << ")" << std::endl;
                    // diagonal of diagonal blocks
                    K.coeffRef(IDi + ii, IDi + ii) += K_e(idi + ii, idi + ii);
                }
                // extradiagonals of diagonal blocks
                K.coeffRef(IDi, IDi + 1) += K_e(idi, idi + 1);
                K.coeffRef(IDi, IDi + 2) += K_e(idi, idi + 2);
                K.coeffRef(IDi + 1, IDi + 2) += K_e(idi + 1, idi + 2);
                K.coeffRef(IDi + 1, IDi) += K_e(idi + 1, idi);
                K.coeffRef(IDi + 2, IDi) += K_e(idi + 2, idi);
                K.coeffRef(IDi + 2, IDi + 1) += K_e(idi + 2, idi + 1);
                for (std::size_t j = i + 1; j < npdofs; j++)
                {
                    std::size_t IDj = pdofs.getID(j) * 3;
                    std::size_t idj = 3 * j;

                    for (std::size_t ii = 0; ii < 3; ii++)
                    {
                        // diagonal of extradiagonal blocks
                        K.coeffRef(IDi + ii, IDj + ii) += K_e(idi + ii, idj + ii);
                        K.coeffRef(IDj + ii, IDi + ii) += K_e(idj + ii, idi + ii);
                    }
                    // extradiagonals of extradiagonal blocks
                    K.coeffRef(IDi, IDj + 1) += K_e(idi, idj + 1);
                    K.coeffRef(IDi, IDj + 2) += K_e(idi, idj + 2);
                    K.coeffRef(IDi + 1, IDj + 2) += K_e(idi + 1, idj + 2);
                    K.coeffRef(IDi + 1, IDj) += K_e(idi + 1, idj);
                    K.coeffRef(IDi + 2, IDj) += K_e(idi + 2, idj);
                    K.coeffRef(IDi + 2, IDj + 1) += K_e(idi + 2, idj + 1);

                    K.coeffRef(IDj, IDi + 1) += K_e(idj, idi + 1);
                    K.coeffRef(IDj, IDi + 2) += K_e(idj, idi + 2);
                    K.coeffRef(IDj + 1, IDi + 2) += K_e(idj + 1, idi + 2);
                    K.coeffRef(IDj + 1, IDi) += K_e(idj + 1, idi);
                    K.coeffRef(IDj + 2, IDi) += K_e(idj + 2, idi);
                    K.coeffRef(IDj + 2, IDi + 1) += K_e(idj + 2, idi + 1);
                }

                if (order == 1)
                {
                    F[IDi] += f[0];
                    F[IDi + 1] += f[1];
                    F[IDi + 2] += f[2];
                }
            }

            for (std::size_t i = 0; i < pdofs.getnumPdofs(); i++)
            {
                F[3 * pdofs.getID(npdofs - pdofs.getnumPdofs() + i)] += f[i];
                F[3 * pdofs.getID(npdofs - pdofs.getnumPdofs() + i) + 1] += f[i + 1];
                F[3 * pdofs.getID(npdofs - pdofs.getnumPdofs() + i) + 2] += f[i + 2];
            }
            */
            end = std::chrono::high_resolution_clock::now();
            // std::cout << "place K_e in K: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
            place_K_e += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        }
        std::cout << "compute_K_e " << compute_K_e << std::endl;
        std::cout << "place_K_e " << place_K_e << std::endl;

        // std::cout<<std::setprecision(16);
        // std::cout<<"F"<<std::endl<<F<<std::endl<<std::endl;
    }
    /*
        void
        enforceHomogeneousDirichletBC(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
                                      const VirtualDofsCollection &DOFS,
                                      const std::vector<std::function<real(real, real, real)>> &constraintFunctions)
        {

            const auto nVDOFS = DOFS.getnumVdofs();
            const auto nEDOFS = DOFS.getnumEdofs();
            const auto nFDOFS = DOFS.getnumFdofs();
            const auto nPDOFS = DOFS.getnumPdofs();
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
                        constrainedDofs.push_back(3 * j);
                        constrainedDofs.push_back(3 * j + 1);
                        constrainedDofs.push_back(3 * j + 2);
                        // add value 0.0 to the constrained dof values vector
                        constrainedDofsValues.push_back(0.0);
                        constrainedDofsValues.push_back(0.0);
                        constrainedDofsValues.push_back(0.0);
                        break;
                    }
                }
                if (unconstrained)
                {
                    // add the dof id to the unconstrained dofs vector
                    unconstrainedDofs.push_back(3 * j);
                    unconstrainedDofs.push_back(3 * j + 1);
                    unconstrainedDofs.push_back(3 * j + 2);
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
                        constrainedDofs.push_back(3 * j);
                        constrainedDofs.push_back(3 * j + 1);
                        constrainedDofs.push_back(3 * j + 2);
                        // add value 0.0 to the constrained dof values vector
                        constrainedDofsValues.push_back(0.0);
                        constrainedDofsValues.push_back(0.0);
                        constrainedDofsValues.push_back(0.0);
                        break;
                    }
                }
                if (unconstrained)
                {
                    // add the dof id to the unconstrained dofs vector
                    unconstrainedDofs.push_back(3 * j);
                    unconstrainedDofs.push_back(3 * j + 1);
                    unconstrainedDofs.push_back(3 * j + 2);
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
                        constrainedDofs.push_back(3 * j);
                        constrainedDofs.push_back(3 * j + 1);
                        constrainedDofs.push_back(3 * j + 2);
                        // add value 0.0 to the constrained dof values vector
                        constrainedDofsValues.push_back(0.0);
                        constrainedDofsValues.push_back(0.0);
                        constrainedDofsValues.push_back(0.0);
                        break;
                    }
                }
                if (unconstrained)
                {
                    // add the dof id to the unconstrained dofs vector
                    unconstrainedDofs.push_back(3 * j);
                    unconstrainedDofs.push_back(3 * j + 1);
                    unconstrainedDofs.push_back(3 * j + 2);
                }
            }
            for (std::size_t p = 0; p < nPDOFS; p++)
            {
                auto j = nVDOFS + nEDOFS + nFDOFS + p;
                // add the dof id to the unconstrained dofs vector
                unconstrainedDofs.push_back(3 * j);
                unconstrainedDofs.push_back(3 * j + 1);
                unconstrainedDofs.push_back(3 * j + 2);
            }
        }
    */
    void
    enforceHomogeneousDirichletBC(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
                                  const VirtualDofsCollection &DOFS,
                                  const std::vector<std::function<real(real, real, real)>> &constraintFunctions)
    {
        const auto nVDOFS = DOFS.getnumVdofs();
        const auto nEDOFS = DOFS.getnumEdofs();
        const auto nFDOFS = DOFS.getnumFdofs();
        const auto nPDOFS = DOFS.getnumPdofs();
        isConstrained.resize(3 * DOFS.getnumDofs(), false);
        std::cout << "size of isConstrained " << isConstrained.size() << std::endl;

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
        for (std::size_t p = 0; p < nPDOFS; p++)
        {
            auto j = nVDOFS + nEDOFS + nFDOFS + p;
            // add the dof id to the unconstrained dofs vector
            /*
            unconstrainedDofs.push_back(3 * j);
            unconstrainedDofs.push_back(3 * j + 1);
            unconstrainedDofs.push_back(3 * j + 2);
            */
        }
        /*
        for (const auto &v : constrainedDofs)
        {
            std::cout << v << std::endl;
        }
        std::cout << "size constr " << constrainedDofs.size() << std::endl;
        std::cout << "size unconstr " << unconstrainedDofs.size() << std::endl;
        std::cout << DOFS.getnumDofs() << std::endl;
        */
    }

    real computeStrainError(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
                            const LocalVirtualDofsCollection &dofs,
                            const VirtualPolyhedronProjections &vp,
                            const std::function<std::array<real, 6>(real, real, real)> &EpsEx_func)
    {
        real error = 0.0;
        unsigned int order = vp.getOrder();
        auto m_kminus1 = Monomial3D::getMonomialsOrdered(order - 1);
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
            // std::cout << C.rows() << " " << C.cols() << " " << U_loc.size() << std::endl;
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
};

#endif // __SOLVER_HPP_