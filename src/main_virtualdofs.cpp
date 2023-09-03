#include "virtualDofs.hpp"
#include <iostream>

int main()
{
    std::string filename = "voro8.geo";
    Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> mesh(filename);

    VirtualDofsCollection DOFS(mesh, 2);
    std::cout << DOFS << std::endl;

    DOFS.print();

    std::cout<<std::endl;

    LocalVirtualDofsCollection dofs(mesh,DOFS);

    return 0;
}