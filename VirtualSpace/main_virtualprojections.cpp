#include "virtualProjections.hpp"
#include <iostream>

int main()
{
    std::string filename = "voro8.geo";
    Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> mesh(filename);

    // std::cout << "here!" << std::endl;
    // std::cout << integrateMonomial(2, mesh.getPolygon(1), Monomial2D(2, 0, 1.0), Point3D(0.5, 0.25, 0.25), sqrt(2) / 2.0, Point3D(0, 0, 1), Point3D(0, 1, 0)) << std::endl;
    // std::cout << integrateMonomial3DRestrictedMonomial2D(Point3D(0.25, 0.25, 0.25), sqrt(3) / 2.0, mesh.getPolygon(1), Monomial3D(0, 0, 2, 1.0), Monomial2D(0, 0, 1.0)) << std::endl;
/*
    Point3D p0(0, 0, 0);
    Point3D p1(3, 0, 0);
    Point3D p2(3, 2, 0);
    Point3D p3(1.5, 4, 0);
    Point3D p4(0, 4, 0);

    Edge3D e1(p0, p1);
    Edge3D e2(p1, p2);
    Edge3D e3(p2, p3);
    Edge3D e4(p3, p4);
    Edge3D e5(p4, p0);

    Point3D p10(2, 0, 0);
    Point3D p20(0, 2, 0);

    Edge3D e10(p0, p10);
    Edge3D e20(p10, p20);
    Edge3D e30(p20, p0);

    // Polygon3D F{e1,e2,e3,e4,e5};
    Polygon3D F{e10, e20, e30};

    std::set<Point3D> vertices;
    for (std::size_t e = 0; e < F.numEdges(); e++)
    {
        // std::cout<<p1[f][e][0]<<" "<<p1[f][e][1]<<std::endl;
        vertices.insert(F[e][0]);
    }
    // Print vertices
    for (const auto &v : vertices)
    {
        std::cout << v << std::endl;
    }
    std::cout << std::endl;
*/

    unsigned int k = 2;
    VirtualDofsCollection DOFS(mesh, k);
    VirtualProjections vp(DOFS, mesh, k);
    //vp.computeFaceProjection(DOFS, F, 1, true);

/*
    real I=0.0;
    I+=integrateMonomial(2,mesh.getPolygon(1),Monomial2D(0,0,2.0),Point3D(0.5,0.25,0.25),sqrt(2.0)/2,Point3D(0,0,1),Point3D(0,1,0));
    I+=integrateMonomial(2,mesh.getPolygon(1),Monomial2D(1,0,0.0),Point3D(0.5,0.25,0.25),sqrt(2.0)/2,Point3D(0,0,1),Point3D(0,1,0));
    I+=integrateMonomial(2,mesh.getPolygon(1),Monomial2D(0,1,0.0),Point3D(0.5,0.25,0.25),sqrt(2.0)/2,Point3D(0,0,1),Point3D(0,1,0));
    I+=integrateMonomial(2,mesh.getPolygon(1),Monomial2D(2,0,-12.0),Point3D(0.5,0.25,0.25),sqrt(2.0)/2,Point3D(0,0,1),Point3D(0,1,0));
    I+=integrateMonomial(2,mesh.getPolygon(1),Monomial2D(1,1,0.0),Point3D(0.5,0.25,0.25),sqrt(2.0)/2,Point3D(0,0,1),Point3D(0,1,0));
    I+=integrateMonomial(2,mesh.getPolygon(1),Monomial2D(0,2,-12.0),Point3D(0.5,0.25,0.25),sqrt(2.0)/2,Point3D(0,0,1),Point3D(0,1,0));
    std::cout<<I<<std::endl;
*/
    return 0;
}