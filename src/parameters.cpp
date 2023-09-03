#include "parameters.hpp"

Parameters::Parameters(const char *parametersFilename)
{
    GetPot datafile(parametersFilename);

    // Get mesh file input string
    mesh_file = datafile("mesh/input_string", "N1.geo");

    // Get values from the filename file or assign defaults
    first_lame = datafile("materials/first_lame", 1.);
    second_lame = datafile("materials/second_lame", 1.);
    order = datafile("VEM/order", 1);
    initialize_face_integrals = datafile("integration/initializeFaceIntegrals", 1);

    youngs_modulus = second_lame * (3.0 * first_lame + 2.0 * second_lame) / (first_lame + second_lame);
    poisson_ratio = first_lame / (2.0 * (first_lame + second_lame));

    // Forcing term
    forcing = [this](real x, real y, real z) -> std::array<real, 3>
    {
        return {-0.1 * std::pow(M_PI, 2) * ((first_lame + second_lame) * cos(M_PI * x) * sin(M_PI * y + M_PI * z) - (first_lame + 4 * second_lame) * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z)),
                -0.1 * std::pow(M_PI, 2) * ((first_lame + second_lame) * cos(M_PI * y) * sin(M_PI * z + M_PI * x) - (first_lame + 4 * second_lame) * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z)),
                -0.1 * std::pow(M_PI, 2) * ((first_lame + second_lame) * cos(M_PI * z) * sin(M_PI * x + M_PI * y) - (first_lame + 4 * second_lame) * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z))};
    };

    // Homogeneous Dirichlet boundary conditions
    HomogeneousDirichletBC.emplace_back([](real x, real y, real z)
                                        { return z; });
    HomogeneousDirichletBC.emplace_back([](real x, real y, real z)
                                        { return x; });
    HomogeneousDirichletBC.emplace_back([](real x, real y, real z)
                                        { return y; });
    HomogeneousDirichletBC.emplace_back([](real x, real y, real z)
                                        { return x - 1.0; });
    HomogeneousDirichletBC.emplace_back([](real x, real y, real z)
                                        { return y - 1.0; });
    HomogeneousDirichletBC.emplace_back([](real x, real y, real z)
                                        { return z - 1.0; });

    // Exact solution
    exact_sol = [](real x, real y, real z) -> std::array<real, 3>
    {
        return {0.1 * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z),
                0.1 * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z),
                0.1 * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z)};
    };

    // Exact strains
    exact_strains = [](real x, real y, real z) -> std::array<real, 6>
    {
        return {0.1 * M_PI * cos(M_PI * x) * sin(M_PI * y) * sin(M_PI * z),
                0.1 * M_PI * sin(M_PI * x) * cos(M_PI * y) * sin(M_PI * z),
                0.1 * M_PI * sin(M_PI * x) * sin(M_PI * y) * cos(M_PI * z),
                0.1 * M_PI * (sin(M_PI * x) * cos(M_PI * y) * sin(M_PI * z) + cos(M_PI * x) * sin(M_PI * y) * sin(M_PI * z)),
                0.1 * M_PI * (sin(M_PI * x) * sin(M_PI * y) * cos(M_PI * z) + sin(M_PI * x) * cos(M_PI * y) * sin(M_PI * z)),
                0.1 * M_PI * (sin(M_PI * x) * sin(M_PI * y) * cos(M_PI * z) + cos(M_PI * x) * sin(M_PI * y) * sin(M_PI * z))};
    };
}

void Parameters::print() const
{
    std::cout << "MATERIAL PARAMETERS" << std::endl;
    std::cout << "Young's modulus E               : " << youngs_modulus << std::endl;
    std::cout << "Poisson ratio  nu               : " << poisson_ratio << std::endl;
    std::cout << "First Lamé parameter lambda     : " << first_lame << std::endl;
    std::cout << "Second Lamé parameter mu        : " << second_lame << std::endl;
    std::cout << std::endl;
    std::cout << "VEM PARAMETERS" << std::endl;
    std::cout << "Order k                         : " << order << std::endl;
}