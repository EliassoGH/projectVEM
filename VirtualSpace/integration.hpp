#ifndef __INTEGRATION_HPP_
#define __INTEGRATION_HPP_

#include "traits.hpp"
#include "mesh.hpp"
#include "monomial.hpp"
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iterator>
#include <utility>
#include <type_traits>

namespace GaussLobatto
{
    using traits::real;

    // Forward declaration of computeGaussLobattoXW
    std::pair<std::vector<real>, std::vector<real>> computeGaussLobattoXW(unsigned int n);

    class GaussLobattoCache
    {
    private:
        // Stores the abscissas and weights for n-points Gauss-Lobatto rule in [-1,1]
        static std::map<unsigned int, std::pair<std::vector<real>, std::vector<real>>> XW;

    public:
        static void initialize(unsigned int n)
        {
            XW.insert(std::make_pair(n, computeGaussLobattoXW(n)));
        }

        static const std::pair<std::vector<real>, std::vector<real>> &getCache(unsigned int n)
        {
            if (XW.find(n) == XW.end()) // Check if the cache for n doesn't exist
            {
                initialize(n); // Initialize the cache for n
            }
            return XW.at(n); // Return the cache for n
        }
    };

    // Initialize static XW
    std::map<unsigned int, std::pair<std::vector<real>, std::vector<real>>> GaussLobattoCache::XW;

    // Function to compute the derivative of the Legendre polynomial of order n
    real legendreDerivative(unsigned int n, real x)
    {
        return (n * x * boost::math::legendre_p(n, x) - n * boost::math::legendre_p(n - 1, x)) / (x * x - 1);
    }

    // Returns the floor((n-4)/2+1) abscissas of n-points Gauss Lobatto inside (0,1) for the standard interval [-1,1]
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

    // Returns the n abscissas and weights of n-points Gauss Lobatto for the standard interval [-1,1]
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

    // Method to compute n-points Gauss-Lobatto points and weights on an edge
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
};

namespace IntegrationMonomial
{
    using traits::real;

    static int count = 1;

    // Integrate a (scaled) 2D or 3D monomial respectively over a non-scaled 2D or 3D domain.
    // To get the integral of a (scaled) monomial over the scaled domain simply multiply the result by h^d.
    // Can also perform integration of a (scaled) monomial in 3D over a 2D domain.
    template <typename DomainType, unsigned int d>
    real integrateMonomial(const unsigned int &N, const DomainType &E, const Monomial<d> &monomial, const Point3D &O = Point3D(), const real &h = 1.0, const Point3D &ex = Point3D(1, 0, 0), const Point3D &ey = Point3D(0, 1, 0))
    {
        if (monomial.getCoefficient() == 0.0)
            return 0.0;
        if (N == 0) // E is a Point
        {
            if constexpr (std::is_same_v<DomainType, Point3D>)
            {
                // std::cout << "N==0 AND DomainType is Point3D" << std::endl;
                // std::cout<<"monomial "<<monomial<<" at "<<E<<" is "<<monomial.evaluate(E)<<std::endl;
                // std::cout<<"evaluation"<<count<<std::endl;
                // count++;
                if (d == 3)
                    return monomial.evaluate((E - O) / h);
                else if (d == 2)
                {
                    Point2D p2D = transformTo2D(E, O, ex, ey, (1.0 / h));
                    // std::cout << E << " " << p2D << "   " << ex << " " << ey << " " << h << std::endl;
                    return monomial.evaluate(p2D);
                }
            }
        }
        else if ((0 < N) && (N < d))
        {
            real I = 0.0;
            real di;
            if constexpr (std::is_same_v<DomainType, Polygon3D>)
            {
                Point3D x0 = (E[0][0] - O) / h;
                for (std::size_t e = 0; e < E.numEdges(); e++)
                {
                    di = ((cross(x0 - ((E[e][0] - O) / h), x0 - ((E[e][1] - O) / h))).norm()) / (E[e].getLength() / h);
                    // std::cout << "N=" << N << ", d=" << d << ", e=" << e << ", di=" << di << std::endl;
                    I += di * integrateMonomial(N - 1, E[e], monomial, O, h);
                }
                for (std::size_t i = 0; i < d; i++)
                {
                    if (monomial.getExponents()[i] > 0)
                    {
                        I += x0[i] * integrateMonomial(N, E, monomial.derivative(i), O, h);
                    }
                }
                I /= (N + monomial.getOrder());
            }
            else if constexpr (std::is_same_v<DomainType, Edge3D>)
            {
                if (d == 3)
                {
                    di = E.getLength() / h;
                    // std::cout<<di<<std::endl;
                    I += di * integrateMonomial(N - 1, E[1], monomial, O, h, ex, ey);
                    Point3D x0 = (E[0] - O) / h;
                    for (std::size_t i = 0; i < d; i++)
                    {
                        if (monomial.getExponents()[i] > 0)
                        {
                            I += x0[i] * integrateMonomial(N, E, monomial.derivative(i), O, h, ex, ey);
                        }
                    }
                }
                else if (d == 2)
                {
                    Point2D x1 = transformTo2D(E[0], O, ex, ey, (1.0 / h));
                    Point2D x2 = transformTo2D(E[1], O, ex, ey, (1.0 / h));
                    Point2D dir = transformTo2D(E.getDirection(), Point3D(), ex, ey);
                    real x01(0.0), x02(0.0);
                    if (std::abs(dir[0]) == 1.0)
                        x02 = x1[1];
                    else if (std::abs(dir[1]) == 1.0)
                        x01 = x1[0];
                    else
                    {
                        if (monomial.getExponents()[0] < monomial.getExponents()[1])
                            x01 = ((x1 - x2)[0] / (x2 - x1)[1] * x2[1] + x2[0]);
                        else
                            x02 = ((x1 - x2)[1] / (x2 - x1)[0] * x2[0] + x2[1]);
                    }
                    Point2D x0(x01, x02);
                    I += (x0 - x1).dot(dir) * monomial.evaluate(x1);
                    // count++;
                    I += (x2 - x0).dot(dir) * monomial.evaluate(x2);
                    // count++;
                    for (std::size_t i = 0; i < d; i++)
                    {
                        if ((monomial.getExponents()[i]) > 0 && (x0[i] != 0.0))
                        {
                            I += x0[i] * integrateMonomial(N, E, monomial.derivative(i), O, h, ex, ey);
                        }
                    }
                }
                I /= (N + monomial.getOrder());
            }
            return I;
        }
        else if (N == d)
        {
            real I = 0.0;
            real bi;
            if constexpr (std::is_same_v<DomainType, Polyhedron<Polygon3D>>)
            {
                for (std::size_t f = 0; f < E.numPolygons(); f++)
                {
                    bi = (E[f].getOutwardNormal().dot((E[f][0][0] - O) / h));
                    // std::cout << "N=" << N << ", d=" << d << ", f=" << f << ", bi=" << bi << std::endl;
                    I += bi * integrateMonomial(N - 1, E[f], monomial, O, h);
                }
                I /= (N + monomial.getOrder());
                I *= std::pow(h, d);
            }
            if constexpr (std::is_same_v<DomainType, Polygon3D>)
            {
                for (std::size_t e = 0; e < E.numEdges(); e++)
                {
                    bi = ((cross(E[e].getDirection(), E.getOutwardNormal())).dot((E[e][0] - O) / h));
                    // std::cout << "N=" << N << ", d=" << d << ", f=" << f << ", bi=" << bi << std::endl;
                    I += bi * integrateMonomial(N - 1, E[e], monomial, O, h, ex, ey);
                }
                I /= (N + monomial.getOrder());
                I *= std::pow(h, d);
            }
            // std::cout<<"functions evaluations: "<<count<<std::endl;
            // count=1;
            return I;
        }
        else
            throw std::logic_error("Integration fault.");
        return 0.0;
    }

    Point3D getPolygonCentroid(const Polygon3D &F)
    {
        Point3D X_F(integrateMonomial(2, F, Monomial3D(1, 0, 0, 1.0)),
                    integrateMonomial(2, F, Monomial3D(0, 1, 0, 1.0)),
                    integrateMonomial(2, F, Monomial3D(0, 0, 1, 1.0)));
        return (X_F / F.getArea());
    }

    // Integrate over a polygon
    real integrateMonomial3DRestrictedMonomial2D(const Point3D &X_P, const real &h_P, const Polygon3D &F, const Monomial3D &m3D, const Monomial2D &m2D)
    {
        Point3D X_F = getPolygonCentroid(F);
        real h_F = F.getDiameter();
        Point3D e_x = F.get_e_x();
        Point3D e_y = F.get_e_y();
        Point3D e_z = F.getOutwardNormal();
        /*
        std::cout << X_F << std::endl
                  << h_F << std::endl
                  << e_x << std::endl
                  << e_y << std::endl
                  << e_z << std::endl;
        */
        auto m3Din2D = toPolynomial2D(m3D, h_P, h_F, X_P, X_F, e_x, e_y, e_z);
        auto poly2D = m3Din2D * m2D;
        real I = 0.0;
        for (const auto &monomialPair : poly2D.getPolynomial())
        {
            auto m = monomialPair.second;
            // std::cout<<m<<std::endl;
            I += integrateMonomial(2, F, m, X_F, h_F, e_x, e_y);
        }
        return I;
    }
}

#endif // __INTEGRATION_HPP_