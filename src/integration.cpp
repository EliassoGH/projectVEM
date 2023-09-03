#include "integration.hpp"

namespace GaussLobatto
{
    // Initialize static XW
    std::map<unsigned int, std::pair<std::vector<real>, std::vector<real>>> GaussLobattoCache::XW;

    void GaussLobattoCache::initialize(unsigned int n)
    {
        XW.insert(std::make_pair(n, computeGaussLobattoXW(n)));
    }

    const std::pair<std::vector<real>, std::vector<real>> &GaussLobattoCache::getCache(unsigned int n)
    {
        if (XW.find(n) == XW.end()) // Check if the cache for n doesn't exist
        {
            initialize(n); // Initialize the cache for n
        }
        return XW.at(n); // Return the cache for n
    }

    std::vector<real> computeOneSideGaussLobattoInternalX(unsigned int n)
    {
        if (n < 4)
            return std::vector<real>();

        // Create a functor representing the function to find roots for
        auto legendreDerivativeFunc = [n](real x)
        { return legendreDerivative(n - 1, x); };

        // 30 bits of precision in the result or machine precision of double
        boost::math::tools::eps_tolerance<double> tolerance(30);
        boost::uintmax_t maxIterations = 100;

        // Define an interval epsilon (legendre is not defined in -1 and 1)
        real interval_epsilon = std::numeric_limits<real>::epsilon() * 10.;
        std::vector<real> result;
        result.reserve(((n - 4) / 2) + 1);
        auto intervals = computeOneSideGaussLobattoInternalX(n - 1);
        std::size_t n_interv = (n - 3) / 2 + 1;
        for (std::size_t i = 0; i < n_interv; i++)
        {
            real a, b;
            if (n > 4)
            {
                if ((n % 2) != 0)
                {
                    if (i == 0)
                    {
                        continue;
                    }
                    else if (i == (n_interv - 1))
                    {
                        a = intervals[i - 1];
                        b = 1.0 - interval_epsilon;
                    }
                    else
                    {
                        a = intervals[i - 1];
                        b = intervals[i];
                    }
                }
                else
                {
                    if (i == 0)
                    {
                        a = 0.0 + interval_epsilon;
                        b = intervals[0];
                    }
                    else if (i == (n_interv - 1))
                    {
                        a = intervals[i - 1];
                        b = 1.0 - interval_epsilon;
                    }
                    else
                    {
                        a = intervals[i - 1];
                        b = intervals[i];
                    }
                }
            }
            else
            {
                a = 0.0 + interval_epsilon;
                b = 1.0 - interval_epsilon;
            }
            std::pair<real, real> abscissa = boost::math::tools::toms748_solve(legendreDerivativeFunc, a, b, tolerance, maxIterations);
            result.push_back((abscissa.first + abscissa.second) / 2.0);
        }
        return result;
    }

    std::pair<std::vector<real>, std::vector<real>> computeGaussLobattoXW(unsigned int n)
    {
        std::vector<real> X;
        std::vector<real> W;
        if (n < 2)
            return std::make_pair(std::vector<real>(), std::vector<real>());
        else
        {
            X.push_back(-1.0);
            W.push_back(2.0 / (n * (n - 1.0)));
            auto input = computeOneSideGaussLobattoInternalX(n);
            std::transform(input.rbegin(), input.rend(), std::back_inserter(X), [](real value)
                           { return -value; });
            if ((n % 2) != 0)
            {
                X.push_back(0.0);
            }
            std::copy(input.begin(), input.end(), std::back_inserter(X));
            std::transform(std::next(X.begin()), X.end(), std::back_inserter(W), [n](real value)
                           { return 2.0 / (n * (n - 1.0) * std::pow((boost::math::legendre_p(n - 1, value)), 2)); });
            X.push_back(1.0);
            W.push_back(2.0 / (n * (n - 1.0)));
        }
        return std::make_pair(X, W);
    }

    std::pair<std::vector<Point3D>, std::vector<real>> computeGaussLobattoPointsWOnEdge(const Edge3D &edge, unsigned int n)
    {
        std::vector<Point3D> points;
        std::vector<real> weights;

        // Get the endpoints of the edge
        const Point3D &startPoint = edge[0];
        const Point3D &endPoint = edge[1];
        const real &length = distance(startPoint, endPoint);

        auto X = GaussLobattoCache::getCache(n).first;
        auto W = GaussLobattoCache::getCache(n).second;

        std::transform(X.begin(), X.end(), std::back_inserter(points), [startPoint, endPoint](real value)
                       { return ((startPoint + endPoint) / 2.0) + value * ((endPoint - startPoint) / 2.0); });
        std::transform(W.begin(), W.end(), std::back_inserter(weights), [length](real value)
                       { return value * (length / 2.0); });

        return std::make_pair(points, weights);
    }
}

namespace IntegrationMonomial
{
    // Initialize static monomials_face_integrals
    std::map<std::size_t, std::vector<real>> MonomialsFaceIntegralsCache::monomials_face_integrals;

    void MonomialsFaceIntegralsCache::initialize(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh, unsigned int order)
    {
        for (const auto &F : mesh.getPolygons())
        {
            initialize(F.second, order);
        }
    }

    void MonomialsFaceIntegralsCache::initialize(const Polygon3D &F, unsigned int order)
    {
        Point3D X_F = getPolygonCentroid(F);
        real h_F = F.getDiameter();
        Point3D e_x = F.get_e_x();
        Point3D e_y = F.get_e_y();
        auto m = Monomial2D::getMonomialsOrdered(order);
        std::vector<real> integrals;
        integrals.reserve(m.size());
        for (std::size_t i = 0; i < m.size(); i++)
        {
            integrals.emplace_back(integrateMonomial(2, F, m[i], X_F, h_F, e_x, e_y));
        }

        auto it = monomials_face_integrals.find(F.getId());
        if (it != monomials_face_integrals.end())
        {
            it->second = integrals;
        }
        else
        {
            monomials_face_integrals.insert(std::make_pair(F.getId(), integrals));
        }
    }

    std::vector<real> &MonomialsFaceIntegralsCache::getCacheMonomials(const Polygon3D &F, unsigned int order)
    {
        auto Id = F.getId();
        auto it = monomials_face_integrals.find(Id);

        if (it == monomials_face_integrals.end()) // Check if the cache for the face doesn't exist
        {
            initialize(F, order); // Initialize the cache
        }
        else if (it->second.size() < ((order + 1) * (order + 2) / 2))
        {
            initialize(F, order);
        }

        return monomials_face_integrals.at(Id); // Return the cache for n
    }

    real &MonomialsFaceIntegralsCache::getCacheMonomial(const Polygon3D &F, const Monomial2D &m)
    {
        auto order = m.getOrder();
        return getCacheMonomials(F, order)[order * (order + 1) / 2 + m.getExponents()[1]];
    }

    Point3D getPolygonCentroid(const Polygon3D &F)
    {
        Point3D X_F(integrateMonomial(2, F, Monomial3D(1, 0, 0, 1.0)),
                    integrateMonomial(2, F, Monomial3D(0, 1, 0, 1.0)),
                    integrateMonomial(2, F, Monomial3D(0, 0, 1, 1.0)));
        return (X_F / F.getArea());
    }

    real integrateMonomial3DRestrictedMonomial2D(const Point3D &X_P, const real &h_P, const Polygon3D &F, const Monomial3D &m3D, const Monomial2D &m2D)
    {
        Point3D X_F = getPolygonCentroid(F);
        real h_F = F.getDiameter();
        Point3D e_x = F.get_e_x();
        Point3D e_y = F.get_e_y();
        Point3D e_z = F.getOutwardNormal();

        auto m3Din2D = toPolynomial2D(m3D, h_P, h_F, X_P, X_F, e_x, e_y, e_z);
        auto poly2D = m3Din2D * m2D;
        auto unit_integrals = MonomialsFaceIntegralsCache::getCacheMonomials(F, poly2D.getOrder());
        real I = 0.0;
        for (const auto &monomialPair : poly2D.getPolynomial())
        {
            auto m = monomialPair.second;
            auto order = m.getOrder();
            if (m.getCoefficient() != 0.0)
                I += unit_integrals[order * (order + 1) / 2 + m.getExponents()[1]] * m.getCoefficient();
        }
        return I;
    }

    real integrateMonomial3DRestrictedPolynomial2D(const Point3D &X_P, const real &h_P, const Polygon3D &F, const Monomial3D &m3D, const Polynomial<2> &p2D, const Point3D &X_F)
    {
        if (X_F == Point3D())
            Point3D X_F = getPolygonCentroid(F);
        real h_F = F.getDiameter();
        Point3D e_x = F.get_e_x();
        Point3D e_y = F.get_e_y();
        Point3D e_z = F.getOutwardNormal();
        real I = 0.0;

        auto m3Din2D = toPolynomial2D(m3D, h_P, h_F, X_P, X_F, e_x, e_y, e_z);
        auto poly2D = m3Din2D * p2D;
        auto unit_integrals = MonomialsFaceIntegralsCache::getCacheMonomials(F, poly2D.getOrder());

        for (const auto &monomialPair : poly2D.getPolynomial())
        {
            auto m = monomialPair.second;
            auto order = m.getOrder();
            if (m.getCoefficient() != 0.0)
                I += unit_integrals[order * (order + 1) / 2 + m.getExponents()[1]] * m.getCoefficient();
        }

        return I;
    }
}

namespace Gauss
{
    std::vector<std::array<real, 4>> gaussPointsTetrahedron(const Point3D &v1, const Point3D &v2, const Point3D &v3, const Point3D &v4, const unsigned int &N_)
    {
        auto v1coords = v1.getCoordinates();
        auto v2coords = v2.getCoordinates();
        auto v3coords = v3.getCoordinates();
        auto v4coords = v4.getCoordinates();
        Eigen::Matrix<real, 4, 3> c;
        c << v1coords[0], v1coords[1], v1coords[2],
            (v2coords[0] - v1coords[0]), (v2coords[1] - v1coords[1]), (v2coords[2] - v1coords[2]),
            (v3coords[0] - v1coords[0]), (v3coords[1] - v1coords[1]), (v3coords[2] - v1coords[2]),
            (v4coords[0] - v1coords[0]), (v4coords[1] - v1coords[1]), (v4coords[2] - v1coords[2]);
        Eigen::Matrix3d subc = c.block<3, 3>(1, 0); // Extract the lower 3x3 submatrix
        real posdet = std::abs(subc.determinant());

        const GaussData &xyzw = data[N_ - 1];

        std::vector<real> W;
        for (unsigned int i = 0; i < xyzw.w.size(); ++i)
        {
            W.emplace_back(xyzw.w[i] * posdet);
        }

        unsigned int NP = std::pow(N_, 3);
        Eigen::MatrixXd xyzmat(NP, 4);
        for (unsigned int i = 0; i < NP; ++i)
        {
            xyzmat.row(i) << 1.0, xyzw.x[i], xyzw.y[i], xyzw.z[i];
        }
        Eigen::MatrixXd XYZ(NP, 3);
        XYZ = xyzmat * c;

        std::vector<std::array<real, 4>> result;
        result.reserve(NP);
        for (Eigen::Index i = 0; i < XYZ.rows(); i++)
        {
            std::array<real, 4> arr = {XYZ(i, 0),
                                       XYZ(i, 1),
                                       XYZ(i, 2),
                                       W[i]};
            result.push_back(arr);
        }

        return result;
    }

    std::vector<std::array<real, 4>> gaussPointsPolyhedron(const Polyhedron<Polygon3D> &P, const unsigned int &N_)
    {
        std::vector<std::array<real, 4>> result;
        Point3D X_P = IntegrationMonomial::getPolyhedronCentroid(P);
        for (std::size_t f = 0; f < P.numPolygons(); f++)
        {
            auto F = P[f];
            Point3D x0 = F[0][0];
            for (std::size_t e = 0; e < (F.numEdges() - 2); e++)
            {
                auto partial = gaussPointsTetrahedron(x0, F[e][1], F[e + 1][1], X_P, N_);
                // std::cout<<x0<<" "<<F[e][1]<<" "<<F[e + 1][1]<<" "<<X_P<<std::endl;
                result.insert(result.end(), partial.begin(), partial.end());
                /*
                std::cout << "partial xyzw" << std::endl;
                for (std::size_t i = 0; i < partial.size(); i++)
                {
                    std::cout << partial[i][0] << " " << partial[i][1] << " " << partial[i][2] << " " << partial[i][3] << " " << std::endl;
                }
                */
            }
        }
        return result;
    }

    real integrateFunctionOverTetrahedron(const std::function<real(real, real, real)> &func,
                                          const Point3D &v1,
                                          const Point3D &v2,
                                          const Point3D &v3,
                                          const Point3D &v4,
                                          unsigned int N_)
    {
        real I = 0.0;
        auto xyzw = gaussPointsTetrahedron(v1, v2, v3, v4, N_);
        for (std::size_t i = 0; i < xyzw.size(); ++i)
        {
            real x = xyzw[i][0];
            real y = xyzw[i][1];
            real z = xyzw[i][2];
            real w = xyzw[i][3];
            I += w * func(x, y, z);
        }
        return I;
    }

    real integrateFunctionOverPolyhedron(const std::function<real(real, real, real)> &func,
                                         const Polyhedron<Polygon3D> &P,
                                         const unsigned int &N_)
    {
        real I = 0.0;
        auto xyzw = gaussPointsPolyhedron(P, N_);
        for (std::size_t i = 0; i < xyzw.size(); ++i)
        {
            real x = xyzw[i][0];
            real y = xyzw[i][1];
            real z = xyzw[i][2];
            real w = xyzw[i][3];
            I += w * func(x, y, z);
        }
        return I;
    }
}