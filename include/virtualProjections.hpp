#ifndef __VIRTUALPROJECTIONS_HPP_
#define __VIRTUALPROJECTIONS_HPP_

#include "mesh.hpp"
#include "parameters.hpp"
#include "monomial.hpp"
#include "integration.hpp"
#include "virtualDofs.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <omp.h>
#include <map>
#include <vector>
#include <chrono>
#include <functional>

using namespace IntegrationMonomial;
using namespace Gauss;

class VirtualFaceProjections
{
private:
    static constexpr double threshold = 1e-12; // entries for matrix G
    std::map<std::size_t, std::vector<Polynomial<2>>> faceProjections;

public:
    /**
     * @brief Constructor
     * 
     * @param dofs VirtualDofsCollection
     * @param mesh 
     * @param order 
     */
    VirtualFaceProjections(const VirtualDofsCollection &dofs, const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh, const unsigned int &order);

    /**
     * @brief Compute the projections of the basis functions corresponding to the dofs defined on the face
     * 
     * @param dofs VirtualDofsCollection
     * @param face face over which the projection is performed
     * @param order 
     * @param checkConsistency if true the code checks the Frobenius norm of the difference between
     *                         the matrix Gf and matrix BfDf, exploiting the identity BfDf=Gf
     * @return std::vector<Polynomial<2>> 
     */
    std::vector<Polynomial<2>>
    computeFaceProjection(const VirtualDofsCollection &dofs, const Polygon3D &face, const unsigned int &order, bool checkConsistency = false);

    /**
     * @brief Get the face projections
     * 
     * @param Id 
     * @return const std::vector<Polynomial<2>>& 
     */
    const std::vector<Polynomial<2>> &getFaceProjection(const std::size_t &Id) const;
};

class VirtualPolyhedronProjections
{
private:
    unsigned int order;
    real youngs_mod;
    real poisson_ratio;
    static constexpr real threshold = 1e-12;                                   // entries for matrix G
    std::map<std::size_t, Eigen::SparseMatrix<real>> G_inv_matrices;           // stores matrices G^-1
    std::map<std::size_t, Eigen::SparseMatrix<real>> polyhedronProjections;    // stores matrices C
    std::map<std::size_t, Eigen::SparseMatrix<real>> elastic_matrices;         // stores matrices E, upper part
    std::map<std::size_t, Eigen::SparseMatrix<real>> deformation_RBM_matrices; // stores matrices Tdr
    std::map<std::size_t, Eigen::VectorXd> forcingProjections;                 // stores polynomials b_hat*V_P

public:
    /**
     * @brief Constructor with the materials parameters and the forcing function
     * 
     * @param E youngs modulus
     * @param nu poisson ratio
     * @param faceProjections
     * @param dofs VirtualDofsCollection
     * @param mesh 
     * @param funcx forcing function x direction
     * @param funcy forcing function y direction
     * @param funcz forcing function z direction
     * @param order 
     */
    VirtualPolyhedronProjections(const real &E, const real &nu, const VirtualFaceProjections &faceProjections, const LocalVirtualDofsCollection &dofs, const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh, const std::function<real(real, real, real)> &funcx, const std::function<real(real, real, real)> &funcy, const std::function<real(real, real, real)> &funcz, const unsigned int &order);

    /**
     * @brief Constructor taking an instance of parameters
     * 
     * @param parameters 
     * @param faceProjections 
     * @param dofs VirtualDofsCollection
     * @param mesh 
     */
    VirtualPolyhedronProjections(const Parameters &parameters, const VirtualFaceProjections &faceProjections, const LocalVirtualDofsCollection &dofs, const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh);

    /**
     * @brief Compute the projections over the polyhedron
     * 
     * @param faceProjections 
     * @param dofs VirtualDofsCollection
     * @param polyhedron 
     * @param funcx forcing function x direction
     * @param funcy forcing function y direction
     * @param funcz forcing function z direction
     * @param order 
     * @return Eigen::SparseMatrix<real> 
     */
    Eigen::SparseMatrix<real> computePolyhedronProjections(const VirtualFaceProjections &faceProjections, const LocalVirtualDofs &dofs, const Polyhedron<Polygon3D> &polyhedron, const std::function<real(real, real, real)> &funcx, const std::function<real(real, real, real)> &funcy, const std::function<real(real, real, real)> &funcz, const unsigned int &order);

    /**
     * @brief Get the matrices G
     * 
     * @return const std::map<std::size_t, Eigen::SparseMatrix<real>>& 
     */
    const std::map<std::size_t, Eigen::SparseMatrix<real>> &getInverseMatricesG() const;

    /**
     * @brief Get the projections C
     * 
     * @return const std::map<std::size_t, Eigen::SparseMatrix<real>>& 
     */
    const std::map<std::size_t, Eigen::SparseMatrix<real>> &getPolyhedronProjections() const;

    /**
     * @brief Get the elastic matrices E
     * 
     * @return const std::map<std::size_t, Eigen::SparseMatrix<real>>& 
     */
    const std::map<std::size_t, Eigen::SparseMatrix<real>> &getElasticMatrices() const;

    /**
     * @brief Get the deformation and rigid body motion matrices Tdr
     * 
     * @return const std::map<std::size_t, Eigen::SparseMatrix<real>>& 
     */
    const std::map<std::size_t, Eigen::SparseMatrix<real>> &getDeformationRBMMatrices() const;

    /**
     * @brief Get the forcing projections
     * 
     * @return const std::map<std::size_t, Eigen::VectorXd>& 
     */
    const std::map<std::size_t, Eigen::VectorXd> &getforcingProjections() const;

    /**
     * @brief Get the order
     * 
     * @return const unsigned int& 
     */
    inline const unsigned int &getOrder() const
    {
        return order;
    }
};

#endif // __VIRTUALPROJECTIONS_HPP_