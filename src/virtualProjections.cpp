#include "virtualProjections.hpp"

VirtualFaceProjections::VirtualFaceProjections(const VirtualDofsCollection &dofs, const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh, const unsigned int &order)
{
    for (const auto &F : mesh.getPolygons())
    {
        faceProjections.insert(std::make_pair(F.second.getId(), computeFaceProjection(dofs, F.second, order)));
    }
}

std::vector<Polynomial<2>>
VirtualFaceProjections::computeFaceProjection(const VirtualDofsCollection &dofs, const Polygon3D &face, const unsigned int &order, bool checkConsistency)
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

    Eigen::SparseMatrix<real> G_F(m.size(), m.size());
    Eigen::MatrixXd B_F(m.size(), nDof);
    B_F.setZero(); // Initialize all elements to zero
    if (order == 1)
    {
        Point2D point_sum(0, 0);
        for (unsigned int i = 0; i < nDof; i++)
        {
            point_sum = point_sum + transformTo2D(face[i][0], X_F, e_x, e_y, 1.0 / h_F);
        }
        G_F.insert(0, 0) = 1.0;
        G_F.insert(0, 1) = point_sum[0] / nDof;
        G_F.insert(0, 2) = point_sum[1] / nDof;

        for (unsigned int j = 0; j < nDof; j++)
        {
            real lengthEplus = face.getEdge(j % nDof).getLength() / h_F;
            real lengthEminus = face.getEdge((j + nDof - 1) % nDof).getLength() / h_F;
            Point2D nEplus = transformTo2D(cross(face.getEdge(j % nDof).getDirection(), face.getOutwardNormal()), Point3D(), e_x, e_y);
            Point2D nEminus = transformTo2D(cross(face.getEdge((j + nDof - 1) % nDof).getDirection(), face.getOutwardNormal()), Point3D(), e_x, e_y);
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
            for (unsigned int ie = 0; ie < (order - 1); ie++)
            {
                auto j = face.numEdges() + (order - 1) * e + ie;
                B_F(0, j) = 0.0;
                auto edgeDofId = dofs.getnumVdofs() + (order - 1) * (edgeId - 1) + ie;
                auto edgeDof = dofs.getDof<EdgeDof>(edgeDofId);
                Point2D xi = transformTo2D(edgeDof->getGaussLobattoPoint(), X_F, e_x, e_y, (1.0 / h_F));
                real weight = edgeDof->getWeight() / h_F;
                for (unsigned int i = 1; i < m.size(); i++)
                {
                    B_F(i, j) = weight * ((m[i].dx()).evaluate(xi) * nE[0] +
                                          (m[i].dy()).evaluate(xi) * nE[1]);
                }
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
        G_F.insert(i, i) = (integrateMonomial(2, face, (m[i].dx() * m[i].dx()), X_F, h_F, e_x, e_y) +
                            integrateMonomial(2, face, (m[i].dy() * m[i].dy()), X_F, h_F, e_x, e_y)) /
                           std::pow(h_F, 2);
        for (std::size_t j = i + 1; j < m.size(); j++)
        {
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
        for (std::size_t j = 0; j < m.size(); j++)
        {
            for (std::size_t i = 0; i < face.numEdges(); i++)
            {
                D_F(i, j) = m[j].evaluate(transformTo2D(face.getEdge(i)[0], X_F, e_x, e_y, (1.0 / h_F)));
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
        Eigen::MatrixXd diff = G_F - (B_F * D_F);
        std::cout << "Froebenius norm = " << diff.norm() << std::endl;
    }

    Eigen::MatrixXd S_F(m.size(), nDof);
    Eigen::SparseLU<Eigen::SparseMatrix<real>> solver;
    solver.analyzePattern(G_F);
    solver.factorize(G_F);
    S_F = solver.solve(B_F);

    std::vector<Polynomial<2>> result;
    auto ordered_monomials = Monomial2D::getMonomialsOrdered(order);
    for (std::size_t j = 0; j < nDof; j++)
    {
        Polynomial<2> poly;
        for (std::size_t i = 0; i < m.size(); i++)
        {
            auto m = ordered_monomials[i];
            m.setCoefficient(S_F(i, j));
            poly.addMonomial(m);
        }
        result.emplace_back(poly);
    }
    return result;
}

const std::vector<Polynomial<2>> &VirtualFaceProjections::getFaceProjection(const std::size_t &Id) const
{
    return faceProjections.at(Id);
}

VirtualPolyhedronProjections::VirtualPolyhedronProjections(const real &E, const real &nu, const VirtualFaceProjections &faceProjections, const LocalVirtualDofsCollection &dofs, const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh, const std::function<real(real, real, real)> &funcx, const std::function<real(real, real, real)> &funcy, const std::function<real(real, real, real)> &funcz, const unsigned int &order)
{
    youngs_mod = E;
    poisson_ratio = nu;
    for (const auto &P : mesh.getPolyhedra())
    {
        polyhedronProjections.insert(std::make_pair(P.second.getId(), computePolyhedronProjections(faceProjections, dofs.getLocalDofs(P.second.getId()), P.second, funcx, funcy, funcz, order)));
        std::cout << P.second.getId() << std::endl;
    }
}

VirtualPolyhedronProjections::VirtualPolyhedronProjections(const Parameters &parameters, const VirtualFaceProjections &faceProjections, const LocalVirtualDofsCollection &dofs, const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh)
{
    order = parameters.getOrder();
    youngs_mod = parameters.getYoungsModulus();
    poisson_ratio = parameters.getPoissonRatio();
    auto forcing = parameters.getForcing();
    std::function<real(real, real, real)> funcx = [forcing](real x, real y, real z) -> real
    {
        return forcing(x, y, z)[0];
    };
    std::function<real(real, real, real)> funcy = [forcing](real x, real y, real z) -> real
    {
        return forcing(x, y, z)[1];
    };
    std::function<real(real, real, real)> funcz = [forcing](real x, real y, real z) -> real
    {
        return forcing(x, y, z)[2];
    };
    for (const auto &P : mesh.getPolyhedra())
    {
        polyhedronProjections.insert(std::make_pair(P.second.getId(), computePolyhedronProjections(faceProjections, dofs.getLocalDofs(P.second.getId()), P.second, funcx, funcy, funcz, order)));
    }
}

Eigen::SparseMatrix<real> VirtualPolyhedronProjections::computePolyhedronProjections(const VirtualFaceProjections &faceProjections, const LocalVirtualDofs &dofs, const Polyhedron<Polygon3D> &polyhedron, const std::function<real(real, real, real)> &funcx, const std::function<real(real, real, real)> &funcy, const std::function<real(real, real, real)> &funcz, const unsigned int &order)
{
    // Materials coefficients
    real a = youngs_mod / ((1.0 + poisson_ratio) * (1.0 - 2 * poisson_ratio)) * (1.0 - poisson_ratio);
    real b = youngs_mod / ((1.0 + poisson_ratio) * (1.0 - 2 * poisson_ratio)) * (poisson_ratio);
    real c = youngs_mod / ((1.0 + poisson_ratio) * (1.0 - 2 * poisson_ratio)) * ((1.0 - 2 * poisson_ratio) / 2.0);

    unsigned int nDof = dofs.getnumDofs();
    Point3D X_P = getPolyhedronCentroid(polyhedron);
    real h_P = polyhedron.getDiameter();
    real V_P = getPolyhedronVolume(polyhedron);

    auto m_kminus1 = Monomial3D::getMonomialsOrdered(order - 1);
    auto m_k = Monomial3D::getMonomialsOrdered(order);
    std::size_t n_kminus2(0), nu_kminus2(0);
    if (order > 1)
    {
        n_kminus2 = Monomial2D::getMonomialsOrdered(order - 2).size();
        nu_kminus2 = Monomial3D::getMonomialsOrdered(order - 2).size();
    }
    // Fill matrix
    Eigen::SparseMatrix<real> G(6 * m_kminus1.size(), 6 * m_kminus1.size());
    Eigen::SparseMatrix<real> E(6 * m_kminus1.size(), 6 * m_kminus1.size());
    Eigen::SparseMatrix<real> Q(3 * nu_kminus2, 3 * nu_kminus2);
    // auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < m_kminus1.size(); i++)
    {
        real entry = integrateMonomial(3, polyhedron, m_kminus1[i] * m_kminus1[i], X_P, h_P);
        for (std::size_t ii = 0; ii < 6; ii++)
        {
            G.insert(6 * i + ii, 6 * i + ii) = entry;
        }
        for (std::size_t ii = 0; ii < 3; ii++)
        {
            E.insert(6 * i + ii, 6 * i + ii) = entry * a;
            E.insert(6 * i + ii + 3, 6 * i + ii + 3) = entry * c;
        }
        E.insert(6 * i, 6 * i + 1) = entry * b;
        E.insert(6 * i, 6 * i + 2) = entry * b;
        E.insert(6 * i + 1, 6 * i + 2) = entry * b;

        if ((order > 2) && (i < nu_kminus2))
        {
            for (std::size_t ii = 0; ii < 3; ii++)
            {
                Q.insert(3 * i + ii, 3 * i + ii) = entry;
            }
        }
        for (std::size_t j = i + 1; j < m_kminus1.size(); j++)
        {
            real entry = integrateMonomial(3, polyhedron, m_kminus1[i] * m_kminus1[j], X_P, h_P);
            if (std::abs(entry) > threshold)
            {
                for (std::size_t ii = 0; ii < 6; ii++)
                {
                    G.insert(6 * i + ii, 6 * j + ii) = entry;
                }
                for (std::size_t ii = 0; ii < 3; ii++)
                {
                    E.insert(6 * i + ii, 6 * j + ii) = entry * a;
                    E.insert(6 * i + ii + 3, 6 * j + ii + 3) = entry * c;
                }
                E.insert(6 * i, 6 * j + 1) = entry * b;
                E.insert(6 * i, 6 * j + 2) = entry * b;
                E.insert(6 * i + 1, 6 * j + 2) = entry * b;
                E.insert(6 * i + 1, 6 * j) = entry * b;
                E.insert(6 * i + 2, 6 * j) = entry * b;
                E.insert(6 * i + 2, 6 * j + 1) = entry * b;

                if ((order > 2) && (j < nu_kminus2))
                {
                    for (std::size_t ii = 0; ii < 3; ii++)
                    {
                        Q.insert(3 * i + ii, 3 * j + ii) = entry;
                    }
                }
            }
        }
    }
    // auto end = std::chrono::high_resolution_clock::now();
    // elapsed_GEQ += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    G.finalize();
    G.makeCompressed();
    E.finalize();
    E.makeCompressed();
    Q.finalize();
    Q.makeCompressed();
    elastic_matrices.insert(std::make_pair(polyhedron.getId(), E));

    // Fill vector f_vec
    // start = std::chrono::high_resolution_clock::now();
    Eigen::VectorXd f_vec(3 * nu_kminus2);
    if (order == 1)
    {
        f_vec.resize(3);
        f_vec[0] = integrateFunctionOverPolyhedron(funcx, polyhedron, 2) / dofs.getnumVdofs();
        f_vec[1] = integrateFunctionOverPolyhedron(funcy, polyhedron, 2) / dofs.getnumVdofs();
        f_vec[2] = integrateFunctionOverPolyhedron(funcz, polyhedron, 2) / dofs.getnumVdofs();

        forcingProjections.insert(std::make_pair(polyhedron.getId(), f_vec));
    }
    else if (order == 2)
    {
        f_vec[0] = integrateFunctionOverPolyhedron(funcx, polyhedron, order + 1);
        f_vec[1] = integrateFunctionOverPolyhedron(funcy, polyhedron, order + 1);
        f_vec[2] = integrateFunctionOverPolyhedron(funcz, polyhedron, order + 1);

        forcingProjections.insert(std::make_pair(polyhedron.getId(), f_vec));
    }
    else // order > 2
    {
        for (std::size_t i = 0; i < nu_kminus2; i++)
        {
            const auto &m = m_k[i];
            auto mfuncx = [&funcx, &m, &X_P, &h_P](real x, real y, real z)
            {
                return funcx(x, y, z) * m.evaluate((Point3D(x, y, z) - X_P) / h_P);
            };
            auto mfuncy = [&funcy, &m, &X_P, &h_P](real x, real y, real z)
            {
                return funcy(x, y, z) * m.evaluate((Point3D(x, y, z) - X_P) / h_P);
            };
            auto mfuncz = [&funcz, &m, &X_P, &h_P](real x, real y, real z)
            {
                return funcz(x, y, z) * m.evaluate((Point3D(x, y, z) - X_P) / h_P);
            };
            f_vec[3 * i] = integrateFunctionOverPolyhedron(mfuncx, polyhedron, order + 1) * V_P;
            f_vec[3 * i + 1] = integrateFunctionOverPolyhedron(mfuncy, polyhedron, order + 1) * V_P;
            f_vec[3 * i + 2] = integrateFunctionOverPolyhedron(mfuncz, polyhedron, order + 1) * V_P;
        }
        Eigen::VectorXd b_hat(3 * nu_kminus2);
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<real>, Eigen::Upper> solverf;
        b_hat = solverf.compute(Q).solve(f_vec);
        forcingProjections.insert(std::make_pair(polyhedron.getId(), b_hat));
    }
    // end = std::chrono::high_resolution_clock::now();
    // elapsed_f += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    // Fill matrix A1 and Tdr
    std::vector<Eigen::Triplet<real>> tripletList;
    tripletList.reserve(6 * m_kminus1.size() * 3 * nDof * polyhedron.numPolygons());

    // start = std::chrono::high_resolution_clock::now();
    Eigen::SparseMatrix<real> A(6 * m_kminus1.size(), 3 * nDof);
    Eigen::SparseMatrix<real> Tdr(3 * nDof, 3 * m_k.size());
    for (std::size_t f = 0; f < polyhedron.numPolygons(); f++)
    {
        auto n = polyhedron.getPolygon(f).getOutwardNormal();
        auto F = (polyhedron.getPolygon(f).getId() > 0) ? polyhedron.getPolygon(f)
                                                        : polyhedron.getPolygon(f).getOtherPolygon();
        auto Fproj = faceProjections.getFaceProjection(F.getId());
        auto X_F = getPolygonCentroid(F);
        real A_F = F.getArea();
        for (std::size_t e = 0; e < F.numEdges(); e++)
        {
            // Vertex contribution
            std::size_t j = dofs.VToLocalId(F[e][0].getId()); // Local Id of the vertex
            for (std::size_t i = nu_kminus2; i < m_kminus1.size(); i++)
            {
                real integral = integrateMonomial3DRestrictedPolynomial2D(X_P, h_P, F, m_kminus1[i], Fproj[e], X_F);
                if (std::abs(integral) > threshold)
                {
                    tripletList.push_back(Eigen::Triplet<real>(6 * i, 3 * j, integral * n[0]));
                    tripletList.push_back(Eigen::Triplet<real>(6 * i + 1, 3 * j + 1, integral * n[1]));
                    tripletList.push_back(Eigen::Triplet<real>(6 * i + 2, 3 * j + 2, integral * n[2]));
                    tripletList.push_back(Eigen::Triplet<real>(6 * i + 3, 3 * j, integral * n[1]));
                    tripletList.push_back(Eigen::Triplet<real>(6 * i + 3, 3 * j + 1, integral * n[0]));
                    tripletList.push_back(Eigen::Triplet<real>(6 * i + 4, 3 * j + 1, integral * n[2]));
                    tripletList.push_back(Eigen::Triplet<real>(6 * i + 4, 3 * j + 2, integral * n[1]));
                    tripletList.push_back(Eigen::Triplet<real>(6 * i + 5, 3 * j, integral * n[2]));
                    tripletList.push_back(Eigen::Triplet<real>(6 * i + 5, 3 * j + 2, integral * n[0]));
                }
            }
        }
        if (order > 1)
        {
            //  Edge contribution
            for (std::size_t e = 0; e < F.numEdges(); e++)
            {
                // Local dof Ids of the edge
                auto EdofsId = dofs.EToLocalId(F.getPositiveEdge(e).getId());
                // Edge dof
                for (std::size_t ee = 0; ee < EdofsId.size(); ee++)
                {
                    std::size_t j = EdofsId[ee];
                    for (std::size_t i = nu_kminus2; i < m_kminus1.size(); i++)
                    {
                        real integral = integrateMonomial3DRestrictedPolynomial2D(X_P, h_P, F, m_kminus1[i], Fproj[F.numEdges() + (order - 1) * e + ee], X_F);
                        if (std::abs(integral) > threshold)
                        {
                            tripletList.push_back(Eigen::Triplet<real>(6 * i, 3 * j, integral * n[0]));
                            tripletList.push_back(Eigen::Triplet<real>(6 * i + 1, 3 * j + 1, integral * n[1]));
                            tripletList.push_back(Eigen::Triplet<real>(6 * i + 2, 3 * j + 2, integral * n[2]));
                            tripletList.push_back(Eigen::Triplet<real>(6 * i + 3, 3 * j, integral * n[1]));
                            tripletList.push_back(Eigen::Triplet<real>(6 * i + 3, 3 * j + 1, integral * n[0]));
                            tripletList.push_back(Eigen::Triplet<real>(6 * i + 4, 3 * j + 1, integral * n[2]));
                            tripletList.push_back(Eigen::Triplet<real>(6 * i + 4, 3 * j + 2, integral * n[1]));
                            tripletList.push_back(Eigen::Triplet<real>(6 * i + 5, 3 * j, integral * n[2]));
                            tripletList.push_back(Eigen::Triplet<real>(6 * i + 5, 3 * j + 2, integral * n[0]));
                        }
                    }
                }
            }

            //  Face contribution
            //  Local dof Ids of the face
            auto FdofsId = dofs.FToLocalId(F.getId());
            // Face dof
            for (std::size_t ff = 0; ff < n_kminus2; ff++)
            {
                std::size_t j = FdofsId[ff];
                for (std::size_t i = 0; i < m_kminus1.size(); i++)
                {
                    real integral = integrateMonomial3DRestrictedPolynomial2D(X_P, h_P, F, m_kminus1[i], Fproj[order * F.numEdges() + ff], X_F);
                    if (std::abs(integral) > threshold)
                    {
                        tripletList.push_back(Eigen::Triplet<real>(6 * i, 3 * j, integral * n[0]));
                        tripletList.push_back(Eigen::Triplet<real>(6 * i + 1, 3 * j + 1, integral * n[1]));
                        tripletList.push_back(Eigen::Triplet<real>(6 * i + 2, 3 * j + 2, integral * n[2]));
                        tripletList.push_back(Eigen::Triplet<real>(6 * i + 3, 3 * j, integral * n[1]));
                        tripletList.push_back(Eigen::Triplet<real>(6 * i + 3, 3 * j + 1, integral * n[0]));
                        tripletList.push_back(Eigen::Triplet<real>(6 * i + 4, 3 * j + 1, integral * n[2]));
                        tripletList.push_back(Eigen::Triplet<real>(6 * i + 4, 3 * j + 2, integral * n[1]));
                        tripletList.push_back(Eigen::Triplet<real>(6 * i + 5, 3 * j, integral * n[2]));
                        tripletList.push_back(Eigen::Triplet<real>(6 * i + 5, 3 * j + 2, integral * n[0]));
                    }
                }

                // Fill Tdr
                for (std::size_t i = 0; i < m_k.size(); i++)
                {
                    auto m2D = dofs.getDof<FaceDof>(j)->getMonomial();
                    real entry = integrateMonomial3DRestrictedMonomial2D(X_P, h_P, F, m_k[i], m2D) / A_F;
                    if (std::abs(entry) > threshold)
                    {
                        Tdr.insert(3 * j, 3 * i) = entry;
                        Tdr.insert(3 * j + 1, 3 * i + 1) = entry;
                        Tdr.insert(3 * j + 2, 3 * i + 2) = entry;
                    }
                }
            }
        }
    }
    // Polyhedron contribution
    for (std::size_t p = 0; p < dofs.getnumPdofs(); p++)
    {
        auto grad = Monomial3D::getGradientsToMonomialsOrdered(order);
        auto j = dofs.PToLocalId(p);
        for (std::size_t i = 0; i < m_kminus1.size(); i++)
        {
            //  If the coefficient of the x component of the gradient is not 0
            //  and the exponents of the monomial correspond to the dof
            if ((grad[i][0].first != 0.0) && (grad[i][0].second == p))
            {
                tripletList.push_back(Eigen::Triplet<real>(6 * i, 3 * j, (-grad[i][0].first * V_P / h_P)));
                tripletList.push_back(Eigen::Triplet<real>(6 * i + 3, 3 * j + 1, (-grad[i][0].first * V_P / h_P)));
                tripletList.push_back(Eigen::Triplet<real>(6 * i + 5, 3 * j + 2, (-grad[i][0].first * V_P / h_P)));
            }
            // If the coefficient of the y component of the gradient is not 0
            // and the exponents of the monomial correspond to the dof
            if ((grad[i][1].first != 0.0) && (grad[i][1].second == p))
            {
                tripletList.push_back(Eigen::Triplet<real>(6 * i + 1, 3 * j + 1, (-grad[i][1].first * V_P / h_P)));
                tripletList.push_back(Eigen::Triplet<real>(6 * i + 3, 3 * j, (-grad[i][1].first * V_P / h_P)));
                tripletList.push_back(Eigen::Triplet<real>(6 * i + 4, 3 * j + 2, (-grad[i][1].first * V_P / h_P)));
            }
            // If the coefficient of the z component of the gradient is not 0
            // and the exponents of the monomial correspond to the dof
            if ((grad[i][2].first != 0.0) && (grad[i][2].second == p))
            {
                tripletList.push_back(Eigen::Triplet<real>(6 * i + 4, 3 * j + 1, (-grad[i][2].first * V_P / h_P)));
                tripletList.push_back(Eigen::Triplet<real>(6 * i + 2, 3 * j + 2, (-grad[i][2].first * V_P / h_P)));
                tripletList.push_back(Eigen::Triplet<real>(6 * i + 5, 3 * j, (-grad[i][2].first * V_P / h_P)));
            }
        }
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    // end = std::chrono::high_resolution_clock::now();
    // elapsed_ATdr += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    A.finalize();
    A.makeCompressed();

    // Fill Tdr Vertices dofs
    std::size_t i0 = 0;
    for (std::size_t i = 0; i < dofs.getnumVdofs(); i++)
    {
        auto VPoint = dofs.getDof<VertexDof>(i)->getVertex();
        for (std::size_t j = 0; j < m_k.size(); j++)
        {
            real entry = m_k[j].evaluate((VPoint - X_P) / h_P);
            Tdr.insert(3 * i, 3 * j) = entry;
            Tdr.insert(3 * i + 1, 3 * j + 1) = entry;
            Tdr.insert(3 * i + 2, 3 * j + 2) = entry;
        }
    }
    // Fill Tdr Edges dofs
    i0 += dofs.getnumVdofs();
    for (std::size_t i = 0; i < dofs.getnumEdofs(); i++)
    {
        auto GLPoint = dofs.getDof<EdgeDof>(i0 + i)->getGaussLobattoPoint();
        for (std::size_t j = 0; j < m_k.size(); j++)
        {
            real entry = m_k[j].evaluate((GLPoint - X_P) / h_P);
            Tdr.insert(3 * (i0 + i), 3 * j) = entry;
            Tdr.insert(3 * (i0 + i) + 1, 3 * j + 1) = entry;
            Tdr.insert(3 * (i0 + i) + 2, 3 * j + 2) = entry;
        }
    }
    // Fill Tdr Polyhedron dofs
    i0 += dofs.getnumEdofs() + dofs.getnumFdofs();
    for (std::size_t i = 0; i < dofs.getnumPdofs(); i++)
    {
        for (std::size_t j = 0; j < m_k.size(); j++)
        {
            real entry = integrateMonomial(3, polyhedron, m_kminus1[i] * m_k[j], X_P, h_P) / V_P; // / std::pow(h_P, 3);
            for (std::size_t ii = 0; ii < 3; ii++)
            {
                Tdr.insert(3 * (i0 + i) + ii, 3 * j + ii) = entry;
            }
        }
    }

    Tdr.finalize();
    Tdr.makeCompressed();
    deformation_RBM_matrices.insert(std::make_pair(polyhedron.getId(), Tdr));

    Eigen::SparseMatrix<real> C(6 * m_kminus1.size(), 3 * nDof);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<real>, Eigen::Upper> solver;
    solver.compute(G);

    Eigen::SparseMatrix<real> I(G.rows(), G.cols());
    I.setIdentity();
    Eigen::SparseMatrix<real> G_inv = solver.solve(I);
    G_inv_matrices.insert(std::make_pair(polyhedron.getId(), G_inv));

    C = solver.solve(A);
    C.prune(threshold, 1);

    return C;
}

const std::map<std::size_t, Eigen::SparseMatrix<real>> &VirtualPolyhedronProjections::getInverseMatricesG() const
{
    return G_inv_matrices;
}

const std::map<std::size_t, Eigen::SparseMatrix<real>> &VirtualPolyhedronProjections::getPolyhedronProjections() const
{
    return polyhedronProjections;
}

const std::map<std::size_t, Eigen::SparseMatrix<real>> &VirtualPolyhedronProjections::getElasticMatrices() const
{
    return elastic_matrices;
}

const std::map<std::size_t, Eigen::SparseMatrix<real>> &VirtualPolyhedronProjections::getDeformationRBMMatrices() const
{
    return deformation_RBM_matrices;
}

const std::map<std::size_t, Eigen::VectorXd> &VirtualPolyhedronProjections::getforcingProjections() const
{
    return forcingProjections;
}