#include "virtualProjections.hpp"
#include <iostream>

int main()
{
    std::string filename = "voro8.geo";
    Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> mesh(filename);

    unsigned int k = 2;
    MonomialsFaceIntegralsCache::initialize(mesh, 2 * k - 1);

    real E = 2.5;
    real nu = 0.25;
    VirtualDofsCollection DOFS(mesh, k);

    std::cout<<"Computing face projections...";
    VirtualFaceProjections vf(DOFS, mesh, k);
    std::cout<<" done."<<std::endl;

    LocalVirtualDofsCollection dofs(mesh, DOFS);

    auto funcx = [](real x, real y, real z)
    {
        return 0.2;
    };
    auto funcy = [](real x, real y, real z)
    {
        return 0.0;
    };
    auto funcz = [](real x, real y, real z)
    {
        return 0.0;
    };

    std::cout<<"Computing polyhedron projections..."<<std::endl;
    VirtualPolyhedronProjections vp(E, nu, vf, dofs, mesh, funcx, funcy, funcz, k);
    std::cout<<"done"<<std::endl;

    return 0;
}