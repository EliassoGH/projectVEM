#include "export_results.hpp"

void plotting::export_results(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh, const std::vector<real> &solution, const char *filename)
{
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
};

void plotting::export_results(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh,
                              const LocalVirtualDofsCollection &dofs,
                              const VirtualPolyhedronProjections &vp,
                              const std::vector<real> &solution,
                              const char *filename)
{
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
    fprintf(fp, "\n");

    // Compute strains and stresses
    std::vector<std::array<real, 6>> strains{mesh.numVertices(), {0., 0., 0., 0., 0., 0.}};
    std::vector<std::array<real, 6>> stresses{mesh.numVertices(), {0., 0., 0., 0., 0., 0.}};
    std::vector<unsigned int> count(mesh.numVertices(), 0);
    const auto m_kminus1 = Monomial3D::getMonomialsOrdered(vp.getOrder() - 1);
    // For every polyhedron
    for (std::size_t p = 0; p < dofs.numLocalDofsCollection(); p++)
    {
        auto P = mesh.getPolyhedron(p);
        Point3D X_P = getPolyhedronCentroid(P);
        real h_P = P.getDiameter();
        const Eigen::SparseMatrix<real> &C = vp.getPolyhedronProjections().at(p);
        const Eigen::SparseMatrix<real> &E = vp.getElasticMatrices().at(p);
        const Eigen::SparseMatrix<real> &G_inv = vp.getInverseMatricesG().at(p);
        auto pdofs = dofs.getLocalDofs(p);
        Eigen::VectorXd U_loc(3 * pdofs.getnumDofs());
        // Find strains coefficients
        for (std::size_t i = 0; i < pdofs.getnumDofs(); i++)
        {
            U_loc[3 * i] = solution[3 * pdofs.getID(i)];
            U_loc[3 * i + 1] = solution[3 * pdofs.getID(i) + 1];
            U_loc[3 * i + 2] = solution[3 * pdofs.getID(i) + 2];
        }
        Eigen::VectorXd EpsApprox_coeff = C * U_loc;
        Eigen::VectorXd SigmaApprox_coeff = G_inv * E * EpsApprox_coeff;

        //  Find vertices
        std::set<Point3D> vertices;
        for (std::size_t f = 0; f < P.numPolygons(); f++)
        {
            const auto &F = P.getPolygon(f);
            for (std::size_t e = 0; e < F.numEdges(); e++)
            {
                vertices.insert(F[e][0]);
            }
        }

        auto findVertexVal = [&m_kminus1, &X_P, &h_P](const Point3D &v, const Eigen::VectorXd &approx_coeff)
        {
            std::array<real, 6> result = {0., 0., 0., 0., 0., 0.};
            // For each of the 6 components of the strain field
            for (std::size_t j = 0; j < 6; j++)
            {
                real approxVal = 0.0;
                // For every monomial of the strain field
                for (std::size_t i = 0; i < m_kminus1.size(); i++)
                {
                    approxVal += approx_coeff[6 * i + j] * m_kminus1[i].evaluate((v - X_P) / h_P);
                }
                result[j] = approxVal;
            };
            return result;
        };

        for (const auto &v : vertices)
        {
            auto &ID = v.getId();
            auto strainsVertex = findVertexVal(v, EpsApprox_coeff);
            auto stressesVertex = findVertexVal(v, SigmaApprox_coeff);
            for (std::size_t i = 0; i < 6; i++)
            {
                strains[ID][i] += strainsVertex[i];
                stresses[ID][i] += stressesVertex[i];
            }
            count[ID]++;
        }
    }

    // Print vertices strains
    fprintf(fp, "FIELD FieldData 2\n");
    fprintf(fp, "Strain 6 %lu double\n", mesh.numVertices());
    for (std::size_t i = 0; i < mesh.numVertices(); i++)
    {
        fprintf(fp, "%g %g %g %g %g %g\n",
                strains[i][0] / count[i],
                strains[i][1] / count[i],
                strains[i][2] / count[i],
                strains[i][3] / count[i],
                strains[i][4] / count[i],
                strains[i][5] / count[i]);
    }
    fprintf(fp, "\n");

    // Print vertices stresses
    fprintf(fp, "Stresses 6 %lu double\n", mesh.numVertices());
    for (std::size_t i = 0; i < mesh.numVertices(); i++)
    {
        fprintf(fp, "%g %g %g %g %g %g\n",
                stresses[i][0] / count[i],
                stresses[i][1] / count[i],
                stresses[i][2] / count[i],
                stresses[i][3] / count[i],
                stresses[i][4] / count[i],
                stresses[i][5] / count[i]);
    }

    fclose(fp);
};