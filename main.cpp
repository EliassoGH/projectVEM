#include "GTMesh/UnstructuredMesh/UnstructuredMesh.h"

int main () {
    // create an instance of a 3D mesh with embedded memory 
    UnstructuredMesh<3, size_t, double> mesh;
    
    // load mesh from vtk file
    //auto meshReaderVTK = mesh.load("meshFile.vtk");
    
    // load mesh from fpma file
    auto meshReaderFPMA = mesh.load("meshFile.fpma");
    
    auto writer = mesh.write("meshFile1.vtk", *meshReaderFPMA, "My first exported mesh using GTMesh");
    return 0;
}