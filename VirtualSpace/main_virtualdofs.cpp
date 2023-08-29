#include "virtualDofs.hpp"
#include <iostream>

int main()
{
    std::string filename = "voro8.geo";
    Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> mesh(filename);

    VirtualDofsCollection DOFS(mesh, 2);
    std::cout << DOFS << std::endl;
/*
    auto a = DOFS.getDof<VertexDof>(0);
    std::cout << a->getId() << " " << a->getVertex() << std::endl;
    std::cout << *a << std::endl;
    std::cout << DOFS.getDof(0)->getId() << std::endl;
*/
    //std::cout << "GL" << DOFS.getDof<EdgeDof>(27)->getGaussLobattoPoint() << std::endl;

    std::cout << DOFS.getnumDofs() << std::endl;
    std::cout << DOFS.getnumVdofs() << std::endl;
    std::cout << DOFS.getnumEdofs() << std::endl;
    std::cout << DOFS.getnumFdofs() << std::endl;
    std::cout << DOFS.getnumPdofs() << std::endl;

    std::cout<<std::endl;
/*
    LocalVirtualDofs dofs(mesh.getPolyhedron(0), DOFS);
    std::cout << dofs.getnumDofs() << std::endl;
    std::cout << dofs.getnumVdofs() << std::endl;
    std::cout << dofs.getnumEdofs() << std::endl;
    std::cout << dofs.getnumFdofs() << std::endl;
    std::cout << dofs.getnumPdofs() << std::endl;
    std::cout << dofs << std::endl;
*/
    //LocalVirtualDofs pdofs(mesh.getPolyhedron(1),DOFS);

    LocalVirtualDofsCollection dofs(mesh,DOFS);
    for(std::size_t i=0;i<mesh.numPolyhedra();i++)
    {
        std::cout<<dofs.getLocalDofs(i)<<std::endl<<std::endl;
    }

    return 0;
}