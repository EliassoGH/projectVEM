#include "problem.hpp"
#include <iostream>

int main()
{
    // Solve the problem contained in parameters.dat and print the L2 strain norm
    Problem problem("parameters.dat", true);

    return 0;
}