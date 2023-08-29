#ifndef __EXPORT_RESULTS_HPP
#define __EXPORT_RESULTS_HPP

#include "mesh.hpp"
//#include <vtkSmartPointer.h>
//#include <vtkPoints.h>
// #include <vtkCellArray.h>
// #include <vtkPolyData.h>
// #include <vtkPolyDataWriter.h>
#include <vector>
#include <fstream>

namespace plotting
{

    void export_results(Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> mesh, std::vector<real> solution, const char *filename = "output.vtk")
    {
        /*
        for (const auto &s : solution)
        {
            std::cout << s << std::endl;
        }
        */

        FILE *fp = fopen(filename, "w");

        fprintf(fp, "# vtk DataFile Version 2.0\n");
        fprintf(fp, filename);
        fprintf(fp, "\nASCII\n");
        fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
        fprintf(fp, "POINTS %lu double\n", mesh.numVertices());

        // Vertices
        for (std::size_t i = 0; i < mesh.numVertices(); i++)
        {
            const auto &V = mesh.getVertex(i);
            fprintf(fp, "%g %g %g\n", V[0], V[1], V[2]);
        }
        fprintf(fp, "\n");

        // Cells
        std::size_t totf(0);
        for (const auto &f : mesh.getPolygons())
        {
            totf += f.second.numEdges();
        }
        fprintf(fp, "CELLS %lu %lu\n", mesh.numVertices() + mesh.numEdges() + mesh.numPolygons(), 2 * mesh.numVertices() + 3 * mesh.numEdges() + mesh.numPolygons() + totf);
        // Points
        for (std::size_t i = 0; i < mesh.numVertices(); i++)
        {
            fprintf(fp, "1 %ld\n", i);
        }
        // Edges
        for (std::size_t i = 1; i < (mesh.numEdges() + 1); i++)
        {
            fprintf(fp, "2 %ld %ld\n", mesh.getEdge(i)[0].getId(), mesh.getEdge(i)[1].getId());
        }
        // Faces
        for (std::size_t i = 1; i < (mesh.numPolygons() + 1); i++)
        {
            const auto &F = mesh.getPolygon(i);
            fprintf(fp, "%ld ", F.numEdges());
            for (std::size_t j = 0; j < (F.numEdges() - 1); j++)
            {
                fprintf(fp, "%ld"
                            " ",
                        F[j][0].getId()); // take first vertex
            }
            fprintf(fp, "%ld"
                        "\n",
                    F[F.numEdges() - 1][0].getId());
        }
        fprintf(fp, "\n");

        // Cell_types
        fprintf(fp, "CELL_TYPES %lu\n", mesh.numVertices() + mesh.numEdges() + mesh.numPolygons());
        // print 1 for all vertices
        for (std::size_t i = 0; i < mesh.numVertices(); i++)
        {
            fprintf(fp, "1\n");
        }
        // print 3 for all edges
        for (std::size_t i = 0; i < mesh.numEdges(); i++)
        {
            fprintf(fp, "3\n");
        }
        // print 7 for all faces (polygons)
        for (std::size_t i = 0; i < mesh.numPolygons(); i++)
        {
            fprintf(fp, "7\n");
        }
        fprintf(fp, "\n");

        // Print vertices displacements
        fprintf(fp, "POINT_DATA %lu\n", mesh.numVertices());
        fprintf(fp, "VECTORS Displacement double\n");
        for (std::size_t i = 0; i < mesh.numVertices(); i++)
        {
            fprintf(fp, "%g %g %g\n", solution[3 * i], solution[3 * i + 1], solution[3 * i + 2]);
        }

        fclose(fp);
        // outputFile.close();

        /*
        // Create VTK points
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for (const Point3D &vertex : displacedVertices)
        {
            points->InsertNextPoint(vertex[0], vertex[1], vertex[2]);
        }

        // Create VTK cells (assuming vertices are connected to form a simple mesh)
        vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
        for (std::size_t i = 0; i < numVertices; ++i)
        {
            vtkIdType cellVertices[1] = {static_cast<vtkIdType>(i)};
            cells->InsertNextCell(1, cellVertices);
        }

        // Create VTK polydata
        vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
        polyData->SetPoints(points);
        polyData->SetVerts(cells);

        // Write VTK file
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetFileName("output.vtk");
        writer->SetInputData(polyData);
        writer->Write();
        */
    };
};

#endif // __EXPORT_RESULTS_HPP