#include "virtualDofs.hpp"
#include <iostream>

int main()
{
    std::string filename = "N2.geo";
    Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> mesh(filename);

    VirtualDofsCollection DOFS(mesh, 4);
    std::cout << DOFS << std::endl;

    auto a = DOFS.getDof<VertexDof>(0);
    std::cout << a->getId() << " " << a->getVertex() << std::endl;
    std::cout << *a << std::endl;
    std::cout << DOFS.getDof(0)->getId() << std::endl;

    std::cout<<"GL"<<DOFS.getDof<EdgeDof>(27)->getGaussLobattoPoint()<<std::endl;

    std::cout<<DOFS.getnumEdofs();

    return 0;
}