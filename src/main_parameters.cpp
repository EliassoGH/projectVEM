#include "parameters.hpp"

using namespace traits;

int main()
{
    Parameters parameters("parameters.dat");
    std::cout << "first Lame: " << parameters.getFirstLame() << std::endl;
    std::cout << "second Lame: " << parameters.getSecondLame() << std::endl;
    std::cout << "Young's modulus: " << parameters.getYoungsModulus() << std::endl;
    std::cout << "Poisson ratio: " << parameters.getPoissonRatio() << std::endl;
    std::cout << "forcing term at (0,0,0): "
              << parameters.getForcing()(0, 0, 0)[0] << ", "
              << parameters.getForcing()(0, 0, 0)[1] << ", "
              << parameters.getForcing()(0, 0, 0)[2] << ", "
              << std::endl;

    return 0;
}