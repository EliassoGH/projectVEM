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
    /**
     * @brief Constructor
     * 
     */
    Problem(const char *parametersFilename = "parameters.dat", bool printError = false);
};

#endif // __PROBLEM_HPP_