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

/**
 * @brief Indexing function for an Eigen::Vector (built-in from Eigen 3.4)
 * 
 * @tparam T First type
 * @tparam T2 Second type (index type)
 * @param full full vector from which extracting the subvector
 * @param ind vector of indices to extract
 * @return Eigen::Matrix<typename T2::Scalar, T::RowsAtCompileTime, T::ColsAtCompileTime, T::Options> 
 */
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
    Eigen::VectorXd U;                          // solution
    std::vector<bool> isConstrained;            // vector storing true if dof is constrained
    std::vector<std::size_t> unconstrainedDofs; // vector storing uncosntrained dofs
    std::vector<std::size_t> constrainedDofs;   // vector storing constrained dofs
    std::vector<real> constrainedDofsValues;    // vector storing constrained dofs

public:
    /**
     * @brief Default constructor
     * 
     */
    Solver() {}

    /**
     * @brief Constructor
     * 
     * @param K 
     * @param F 
     */
    Solver(const Eigen::SparseMatrix<real> &K, const Eigen::VectorXd &F) : K(K), F(F) {}

    /**
     * @brief Solve the system KU=F
     * 
     */
    void solve();

    /**
     * @brief Getter for F
     * 
     * @return const Eigen::VectorXd& 
     */
    const Eigen::VectorXd &getRightHandSide() const;

    /**
     * @brief Getter for U as Eigen::Vector
     * 
     * @return const Eigen::VectorXd& 
     */
    const Eigen::VectorXd &getSolutionDisplacementsEig() const;

    /**
     * @brief Getter for U as std::vector
     * 
     * @return std::vector<real> 
     */
    std::vector<real> getSolutionDisplacements() const;

    real computeH1error(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
                        const std::function<std::array<real, 3>(real, real, real)> &Uex_func);
};

class SolverVEM : public Solver
{
private:
    static constexpr real tolerance = 1e-12; // for Dirichlet boundary conditions

public:
    /**
     * @brief Constructor, assembles K and F
     * 
     * @param parameters 
     * @param mesh 
     * @param DOFS 
     * @param dofs 
     * @param vp 
     */
    SolverVEM(const Parameters &parameters,
              const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
              const VirtualDofsCollection &DOFS,
              const LocalVirtualDofsCollection &dofs,
              const VirtualPolyhedronProjections &vp);

    /**
     * @brief Enforce homogeneous Dirichlet boundary conditions
     * 
     * @param mesh 
     * @param DOFS 
     * @param constraintFunctions 
     */
    void enforceHomogeneousDirichletBC(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
                                       const VirtualDofsCollection &DOFS,
                                       const std::vector<std::function<real(real, real, real)>> &constraintFunctions);

    /**
     * @brief Compute the error in the L2 strain norm
     * 
     * @param mesh 
     * @param dofs 
     * @param vp 
     * @param EpsEx_func 
     * @return real 
     */
    real computeStrainError(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
                            const LocalVirtualDofsCollection &dofs,
                            const VirtualPolyhedronProjections &vp,
                            const std::function<std::array<real, 6>(real, real, real)> &EpsEx_func);
};

#endif // __SOLVER_HPP_