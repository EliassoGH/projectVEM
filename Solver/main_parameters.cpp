#include "parameters.hpp"

using namespace traits;

int main()
{
    Parameters parameters("parameters.dat");
    std::cout << parameters.getFirstLame() << std::endl;
    std::cout << parameters.getSecondLame() << std::endl;
    std::cout << parameters.getYoungsModulus() << std::endl;
    std::cout << parameters.getPoissonRatio() << std::endl;
    std::cout << parameters.getForcing()(0, 0, 0)[0] << std::endl;

    // const std::function<real(real, real, real)> funcx = parameters.getForcing();
    auto forcing = parameters.getForcing();
    std::function<real(real, real, real)> forcingx = [forcing](real x, real y, real z) -> real
    {
        return forcing(x, y, z)[0];
    };

    return 0;
}