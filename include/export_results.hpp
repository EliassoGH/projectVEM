#ifndef __EXPORT_RESULTS_HPP
#define __EXPORT_RESULTS_HPP

#include "mesh.hpp"
#include "virtualProjections.hpp"
#include <vector>
#include <fstream>

namespace plotting
{
    /**
     * @brief Function to write the output file in the .vtk format
     * 
     * @param mesh
     * @param solution 
     * @param filename 
     */
    void export_results(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
                        const std::vector<real> &solution,
                        const char *filename = "output.vtk");

    
    void export_results(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
                        const LocalVirtualDofsCollection &dofs,
                        const VirtualPolyhedronProjections &vp,
                        const std::vector<real> &solution,
                        const char *filename = "output.vtk");
};

#endif // __EXPORT_RESULTS_HPP