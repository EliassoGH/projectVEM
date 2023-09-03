#ifndef __PARAMETERS_HPP_
#define __PARAMETERS_HPP_

#include "traits.hpp"
#include "GetPot"
#include <functional>
#include <vector>
#include <string>

using namespace traits;

class Parameters
{
private:
    std::string mesh_file;
    unsigned int order;
    bool initialize_face_integrals;
    real youngs_modulus;
    real poisson_ratio;
    real first_lame;
    real second_lame;
    std::function<std::array<real, 3>(real, real, real)> forcing;
    std::vector<std::function<real(real, real, real)>> HomogeneousDirichletBC;
    std::function<std::array<real, 3>(real, real, real)> exact_sol;
    std::function<std::array<real, 6>(real, real, real)> exact_strains;

public:
    /**
     * @brief Constructor
     * 
     * @param parametersFilename 
     */
    Parameters(const char *parametersFilename);

    /**
     * @brief Get the mesh string name
     * 
     * @return const std::string& 
     */
    inline const std::string &getInputMesh() const
    {
        return mesh_file;
    }

    /**
     * @brief Get the homogeneous dirichlet boundary conditions
     * 
     * @return const std::vector<std::function<real(real, real, real)>>& 
     */
    inline const std::vector<std::function<real(real, real, real)>> &getHomogeneousDirichletBC() const
    {
        return HomogeneousDirichletBC;
    }

    /**
     * @brief Get the exact displacement solution
     * 
     * @return const std::function<std::array<real, 3>(real, real, real)>& 
     */
    inline const std::function<std::array<real, 3>(real, real, real)> &getExactSolution() const
    {
        return exact_sol;
    }

    /**
     * @brief Get the exact strains
     * 
     * @return const std::function<std::array<real, 6>(real, real, real)>& 
     */
    inline const std::function<std::array<real, 6>(real, real, real)> &getExactStrains() const
    {
        return exact_strains;
    }

    /**
     * @brief Get the forcing term
     * 
     * @return const std::function<std::array<real, 3>(real, real, real)>& 
     */
    inline const std::function<std::array<real, 3>(real, real, real)> &getForcing() const
    {
        return forcing;
    }

    /**
     * @brief Get the order
     * 
     * @return const unsigned int& 
     */
    inline const unsigned int &getOrder() const
    {
        return order;
    }

    /**
     * @brief Get the bool returning true if face integrals are required to be initialized
     * 
     * @return true 
     * @return false 
     */
    inline const bool &getInitializeFaceIntegrals() const
    {
        return initialize_face_integrals;
    }

    /**
     * @brief Get Young's modulus
     * 
     * @return const real& 
     */
    inline const real &getYoungsModulus() const
    {
        return youngs_modulus;
    }

    /**
     * @brief Get Poisson's ratio
     * 
     * @return const real& 
     */
    inline const real &getPoissonRatio() const
    {
        return poisson_ratio;
    }

    /**
     * @brief Get first Lame's parameter
     * 
     * @return const real& 
     */
    inline const real &getFirstLame() const
    {
        return first_lame;
    }

    /**
     * @brief Get second Lame's parameter
     * 
     * @return const real& 
     */
    inline const real &getSecondLame() const
    {
        return second_lame;
    }

    /**
     * @brief Print parameter information
     * 
     */
    void print() const;
};

#endif // __PARAMETERS_HPP_