#ifndef __VIRTUALPROJECTIONS_HPP_
#define __VIRTUALPROJECTIONS_HPP_

#include "mesh.hpp"
#include "monomial.hpp"
#include "integration.hpp"
#include "virtualDofs.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <map>
#include <vector>

using namespace IntegrationMonomial;

class VirtualProjections
{
private:
    std::map<std::size_t, std::vector<std::vector<Monomial2D>>> faceProjections;
    static constexpr double threshold = 1e-12;

public:
    VirtualProjections(const VirtualDofsCollection &dofs, const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh, const unsigned int &order)
    {
        
        for (const auto &F : mesh.getPolygons())
        {
            computeFaceProjection(dofs, F.second, order, true);
            //std::cout<<F.second<<std::endl;
        }
        
        //computeFaceProjection(dofs, mesh.getPolygon(1), order, true);
    }

    void computeFaceProjection(const VirtualDofsCollection &dofs, const Polygon3D &face, const unsigned int &order, bool checkConsistency = false)
    {
        unsigned int nDof = order * face.numEdges();
        if (order > 1)
            nDof += Monomial2D::getMonomialsOrdered(order - 2).size();
        auto m = Monomial2D::getMonomialsOrdered(order);
        Point3D X_F = getPolygonCentroid(face);
        real h_F = face.getDiameter();
        Point3D e_x = face.get_e_x();
        Point3D e_y = face.get_e_y();
        real A_F = face.getArea();
        /*
        std::cout << X_F << std::endl
                  << h_F << std::endl
                  << e_x << std::endl
                  << e_y << std::endl
                  << A_F << std::endl;
        */

        Eigen::SparseMatrix<real> G_F(m.size(), m.size());
        Eigen::MatrixXd B_F(m.size(), nDof);
        B_F.setZero(); // Initialize all elements to zero
        if (order == 1)
        {
            Point2D point_sum(0, 0);
            for (unsigned int i = 0; i < nDof; i++)
            {
                point_sum = point_sum + transformTo2D(face[i][0], X_F, e_x, e_y, 1.0 / h_F);
                // point_sum = point_sum + transformTo2D(face[i][0], X_F, e_x, e_y);
            }
            G_F.insert(0, 0) = 1.0;
            G_F.insert(0, 1) = point_sum[0] / nDof;
            G_F.insert(0, 2) = point_sum[1] / nDof;

            for (unsigned int j = 0; j < nDof; j++)
            {
                real lengthEplus = face.getEdge(j % nDof).getLength() / h_F;
                real lengthEminus = face.getEdge((j + nDof - 1) % nDof).getLength() / h_F;
                // real lengthEplus = face.getEdge(j % nDof).getLength();
                // real lengthEminus = face.getEdge((j + nDof - 1) % nDof).getLength();
                Point2D nEplus = transformTo2D(cross(face.getEdge(j % nDof).getDirection(), face.getOutwardNormal()), Point3D(), e_x, e_y);
                Point2D nEminus = transformTo2D(cross(face.getEdge((j + nDof - 1) % nDof).getDirection(), face.getOutwardNormal()), Point3D(), e_x, e_y);
                // std::cout << nEplus << " " << nEminus << std::endl;
                B_F(0, j) = 1.0 / nDof;
                B_F(1, j) = 0.5 * (lengthEplus * nEplus[0] + lengthEminus * nEminus[0]);
                B_F(2, j) = 0.5 * (lengthEplus * nEplus[1] + lengthEminus * nEminus[1]);
            }
        }
        else
        {
            for (unsigned int j = 0; j < m.size(); j++)
            {
                real entry = integrateMonomial(2, face, m[j], X_F, h_F, e_x, e_y) / A_F;
                if (std::abs(entry) > threshold)
                    G_F.insert(0, j) = entry;
            }
            for (unsigned int j = 0; j < face.numEdges(); j++)
            {
                B_F(0, j) = 0.0;
                real lengthEplus = face.getEdge(j % face.numEdges()).getLength() / h_F;
                real lengthEminus = face.getEdge((j + face.numEdges() - 1) % face.numEdges()).getLength() / h_F;
                Point2D nEplus = transformTo2D(cross(face.getEdge(j % face.numEdges()).getDirection(), face.getOutwardNormal()), Point3D(), e_x, e_y);
                Point2D nEminus = transformTo2D(cross(face.getEdge((j + face.numEdges() - 1) % face.numEdges()).getDirection(), face.getOutwardNormal()), Point3D(), e_x, e_y);
                Point2D xi = transformTo2D(face.getEdge(j)[0], X_F, e_x, e_y, (1.0 / h_F));
                for (unsigned int i = 1; i < m.size(); i++)
                {
                    B_F(i, j) = (1.0 / (order * (order + 1.0))) * ((m[i].dx()).evaluate(xi) * (lengthEplus * nEplus[0] + lengthEminus * nEminus[0]) +
                                                                   (m[i].dy()).evaluate(xi) * (lengthEplus * nEplus[1] + lengthEminus * nEminus[1]));
                }
            }
            for (unsigned int e = 0; e < face.numEdges(); e++)
            {
                auto edge = face.getPositiveEdge(e);
                auto edgeId = edge.getId();
                Point2D nE = transformTo2D(cross(face.getEdge(e).getDirection(), face.getOutwardNormal()), Point3D(), e_x, e_y);
                // std::cout << "edge= " << face.getEdge(e) << ", positive edge= " << edge << std::endl;
                for (unsigned int ie = 0; ie < (order - 1); ie++)
                {
                    auto j = face.numEdges() + (order - 1) * e + ie;
                    B_F(0, j) = 0.0;
                    auto edgeDofId = dofs.getnumVdofs() + (order - 1) * (edgeId - 1) + ie;
                    // std::cout << "index= " << edgeDofId << std::endl;
                    auto edgeDof = dofs.getDof<EdgeDof>(edgeDofId);
                    Point2D xi = transformTo2D(edgeDof->getGaussLobattoPoint(), X_F, e_x, e_y, (1.0 / h_F));
                    real weight = edgeDof->getWeight() / h_F;
                    for (unsigned int i = 1; i < m.size(); i++)
                    {
                        B_F(i, j) = weight * ((m[i].dx()).evaluate(xi) * nE[0] +
                                              (m[i].dy()).evaluate(xi) * nE[1]);
                    }
                    // std::cout << "point: " << xi << "n: " << nE << std::endl;
                }
            }
            B_F(0, order * face.numEdges()) = 1.0;
            auto lapl = Monomial2D::getLaplaciansToMonomialsOrdered(order);
            for (unsigned int j = (order * face.numEdges()); j < nDof; j++)
            {
                auto jf = j - (order * face.numEdges());
                for (unsigned int i = 1; i < m.size(); i++)
                {
                    if ((lapl[i].first.first != 0.0) && (lapl[i].first.second == jf))
                        B_F(i, j) = -(lapl[i].first.first) * A_F / std::pow(h_F, 2);
                    if ((lapl[i].second.first != 0.0) && (lapl[i].second.second == jf))
                        B_F(i, j) = -(lapl[i].second.first) * A_F / std::pow(h_F, 2);
                }
            }
        }
        for (std::size_t i = 1; i < m.size(); i++)
        {
            // std::cout << "(" << i << ", " << i << ")" << m[i].dx() * m[i].dx() <<" + "<< m[i].dy() * m[i].dy() << std::endl;
            G_F.insert(i, i) = (integrateMonomial(2, face, (m[i].dx() * m[i].dx()), X_F, h_F, e_x, e_y) +
                                integrateMonomial(2, face, (m[i].dy() * m[i].dy()), X_F, h_F, e_x, e_y)) /
                               std::pow(h_F, 2);
            for (std::size_t j = i + 1; j < m.size(); j++)
            {
                // std::cout << "(" << i << ", " << j << ")" << m[i].dx() * m[j].dx() <<" + "<< m[i].dy() * m[j].dy() << std::endl;
                real entry = (integrateMonomial(2, face, (m[i].dx() * m[j].dx()), X_F, h_F, e_x, e_y) +
                              integrateMonomial(2, face, (m[i].dy() * m[j].dy()), X_F, h_F, e_x, e_y)) /
                             std::pow(h_F, 2);
                if (std::abs(entry) > threshold)
                {
                    G_F.insert(i, j) = entry;
                    G_F.insert(j, i) = entry;
                }
            }
        }
        G_F.finalize();
        G_F.makeCompressed();
        if (checkConsistency == true)
        {
            Eigen::MatrixXd D_F(nDof, m.size());
            // D_F.setZero(); // Initialize all elements to zero
            for (std::size_t j = 0; j < m.size(); j++)
            {
                for (std::size_t i = 0; i < face.numEdges(); i++)
                {
                    D_F(i, j) = m[j].evaluate(transformTo2D(face.getEdge(i)[0], X_F, e_x, e_y, (1.0 / h_F)));
                    // D_F(i, j) = m[j].evaluate(transformTo2D(face.getEdge(i)[0], X_F, e_x, e_y));
                }
                for (unsigned int e = 0; e < face.numEdges(); e++)
                {
                    for (unsigned int ie = 0; ie < (order - 1); ie++)
                    {
                        auto edge = face.getPositiveEdge(e);
                        auto edgeId = edge.getId();
                        auto i = face.numEdges() + (order - 1) * e + ie;
                        auto edgeDofId = dofs.getnumVdofs() + (order - 1) * (edgeId - 1) + ie;
                        auto edgeDof = dofs.getDof<EdgeDof>(edgeDofId);
                        Point2D xi = transformTo2D(edgeDof->getGaussLobattoPoint(), X_F, e_x, e_y, (1.0 / h_F));
                        D_F(i, j) = m[j].evaluate(xi);
                    }
                }
                for (unsigned int i = (order * face.numEdges()); i < nDof; i++)
                {
                    auto f = i - (order * face.numEdges());
                    D_F(i, j) = integrateMonomial(2, face, (m[j] * m[f]), X_F, h_F, e_x, e_y) / A_F;
                }
            }
/*
            Eigen::MatrixXd denseMatrix = G_F;
            std::cout << denseMatrix << std::endl
                      << std::endl;
            std::cout << B_F << std::endl
                      << std::endl;
            std::cout << D_F << std::endl
                      << std::endl;
            std::cout << (B_F * D_F) << std::endl
                      << std::endl;
            std::cout << denseMatrix << std::endl
                      << std::endl;
*/
            Eigen::MatrixXd diff = G_F - (B_F * D_F);
/*
            for (std::size_t i = 0; i < diff.rows(); i++)
            {
                for (std::size_t j = 0; j < diff.cols(); j++)
                {
                    if (std::abs(diff(i, j)) > 1e-9)
                        std::cout << i << " " << j << std::endl;
                }
            }
*/
            std::cout << "Froebenius norm = " << diff.norm() << std::endl;
        }

        Eigen::MatrixXd S_F(m.size(), nDof);
        Eigen::SparseLU<Eigen::SparseMatrix<real>> solver;
        solver.analyzePattern(G_F);
        solver.factorize(G_F);
        S_F = solver.solve(B_F);
        // S_F = G_F.colPivHouseholderQr().solve(B_F);
        // std::cout << S_F << std::endl;
    }
};

#endif // __VIRTUALPROJECTIONS_HPP_