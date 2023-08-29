#ifndef __INTEGRATION_HPP_
#define __INTEGRATION_HPP_

#include "traits.hpp"
#include "mesh.hpp"
#include "monomial.hpp"
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iterator>
#include <utility>
#include <type_traits>
#include <chrono>
#include <map>
#include <unordered_map>

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

    // Forward declarations
    template <typename DomainType, unsigned int d>
    real integrateMonomial(const unsigned int &N, const DomainType &E, const Monomial<d> &monomial, const Point3D &O = Point3D(), const real &h = 1.0, const Point3D &ex = Point3D(1, 0, 0), const Point3D &ey = Point3D(0, 1, 0));
    Point3D getPolygonCentroid(const Polygon3D &F);

    class MonomialsFaceIntegralsCache
    {
    private:
        // Stores the integrals of the monomials for the faces
        static std::map<std::size_t, std::vector<real>> monomials_face_integrals;
        // static std::vector<std::vector<real>> monomials_face_integrals;

    public:
        static void initialize(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh, unsigned int order)
        {
            for (const auto &F : mesh.getPolygons())
            {
                initialize(F.second, order);
            }
        }
        static void initialize(const Polygon3D &F, unsigned int order)
        {
            // monomials_face_integrals.reserve(231);
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
            // monomials_face_integrals[F.getId()-1]=integrals;

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

        static std::vector<real> &getCacheMonomials(const Polygon3D &F, unsigned int order)
        {
            auto Id = F.getId();
            auto it = monomials_face_integrals.find(Id);

            if (it == monomials_face_integrals.end()) // Check if the cache for the face doesn't exist
            {
                std::cout << "initializing for an order" << std::endl;
                initialize(F, order); // Initialize the cache
            }
            else if (it->second.size() < ((order + 1) * (order + 2) / 2))
            {
                std::cout << "increasing order" << std::endl;
                initialize(F, order);
            }

            // std::cout<<", calling monomial n. "<<((order) * (order + 1) / 2) + m.getExponents()[1]<<std::endl;
            return monomials_face_integrals.at(Id); // Return the cache for n

            // return monomials_face_integrals[Id - 1];
        }

        static real &getCacheMonomial(const Polygon3D &F, const Monomial2D &m)
        {
            auto order = m.getOrder();
            return getCacheMonomials(F, order)[order * (order + 1) / 2 + m.getExponents()[1]];
        }
    };

    // Initialize static monomials_face_integrals
    std::map<std::size_t, std::vector<real>> MonomialsFaceIntegralsCache::monomials_face_integrals;
    // std::vector<std::vector<real>> MonomialsFaceIntegralsCache::monomials_face_integrals={};

    static int count = 1;
    static int64_t elapsed_t = 0;

    // Integrate a (scaled) 2D or 3D monomial respectively over a non-scaled 2D or 3D domain.
    // To get the integral of a (scaled) monomial over the scaled domain simply multiply the result by h^d.
    // Can also perform integration of a (scaled) monomial in 3D over a 2D domain.
    template <typename DomainType, unsigned int d>
    real integrateMonomial(const unsigned int &N, const DomainType &E, const Monomial<d> &monomial, const Point3D &O, const real &h, const Point3D &ex, const Point3D &ey)
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
                /*
                Point3D x0 = (E[0][0] - O) / h;
                for (std::size_t e = 1; e < (E.numEdges() - 1); e++)
                {
                    di = (((E[e][0] - O) / h) - x0).dot(cross(E[e][1] - E[e][0], E.getOutwardNormal()).normalize());
                    // di = ((cross(x0 - ((E[e][0] - O) / h), x0 - ((E[e][1] - O) / h))).norm()) / (E[e].getLength() / h);
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
                */
                Point3D n = E.getOutwardNormal();
                real x01(0.0), x02(0.0), x03(0.0);
                auto exponents = monomial.getExponents();
                // std::cout<<n<<std::endl;
                if (std::abs(n[0]) > 0.99) // x01!=0
                {
                    // std::cout<<"x"<<std::endl;
                    x01 = n.dot((E[0][0] - O) / h) / n[0];
                }
                else if (std::abs(n[1]) > 0.99) // x02!=0
                {
                    // std::cout<<"y"<<std::endl;
                    x02 = n.dot((E[0][0] - O) / h) / n[1];
                }
                else if (std::abs(n[2]) > 0.99) // x03!=0
                {
                    // std::cout<<"z"<<std::endl;
                    x03 = n.dot((E[0][0] - O) / h) / n[2];
                }
                else if (std::abs(n[0]) < 0.1)
                {
                    if (exponents[1] < exponents[2]) // x02!=0
                    {
                        x02 = n.dot((E[0][0] - O) / h) / n[1];
                    }
                    else // x03!=0
                    {
                        x03 = n.dot((E[0][0] - O) / h) / n[2];
                    }
                }
                else if (std::abs(n[1]) < 0.1)
                {
                    if (exponents[0] < exponents[2]) // x01!=0
                    {
                        x01 = n.dot((E[0][0] - O) / h) / n[0];
                    }
                    else // x03!=0
                    {
                        x03 = n.dot((E[0][0] - O) / h) / n[2];
                    }
                }
                else if (std::abs(n[2]) < 0.1)
                {
                    if (exponents[0] < exponents[1]) // x01!=0
                    {
                        x01 = n.dot((E[0][0] - O) / h) / n[0];
                    }
                    else // x02!=0
                    {
                        x02 = n.dot((E[0][0] - O) / h) / n[1];
                    }
                }
                else
                {
                    if (exponents[0] > exponents[1])
                    {
                        if (exponents[1] > exponents[2]) // x03!=0
                        {
                            x03 = n.dot((E[0][0] - O) / h) / n[2];
                        }
                        else // x02!=0
                        {
                            x02 = n.dot((E[0][0] - O) / h) / n[1];
                        }
                    }
                    else
                    {
                        if (exponents[0] < exponents[2]) // x01!=0
                        {
                            x01 = n.dot((E[0][0] - O) / h) / n[0];
                        }
                        else // x03!=0
                        {
                            x03 = n.dot((E[0][0] - O) / h) / n[2];
                        }
                    }
                }
                Point3D x0(x01, x02, x03);
                // std::cout<<x0<<std::endl;
                for (std::size_t e = 0; e < E.numEdges(); e++)
                {
                    di = (((E[e][0] - O) / h) - x0).dot(cross(E[e][1] - E[e][0], n).normalize());
                    // di = ((cross(x0 - ((E[e][0] - O) / h), x0 - ((E[e][1] - O) / h))).norm()) / (E[e].getLength() / h);
                    // std::cout << "N=" << N << ", d=" << d << ", e=" << e << ", di=" << di << std::endl;
                    I += di * integrateMonomial(N - 1, E[e], monomial, O, h);
                }
                for (std::size_t i = 0; i < d; i++)
                {
                    if ((exponents[i] > 0) && (x0[i] != 0.0))
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
                    /*
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
                    */
                    Point3D x1 = (E[0] - O) / h;
                    Point3D x2 = (E[1] - O) / h;
                    Point3D dir = E.getDirection();
                    // std::cout << E[0] << " " << E[1] << " " << O << " " << h << std::endl;
                    // std::cout << x1 << " " << x2 << " " << dir << std::endl;
                    real x01(0.0), x02(0.0), x03(0.0);
                    auto exponents = monomial.getExponents();
                    // std::cout<<"dir: "<<dir<<std::endl;
                    // std::cout<<"abs direction is 1?: "<<(std::abs(dir[0])>0.95)<<" "<<(std::abs(dir[1])>0.95)<<std::endl;
                    if (std::abs(dir[0]) > 0.99) // x01 = 0
                    {
                        // std::cout << "parallel to x" << std::endl;
                        x02 = ((x1 - x2)[1] / (x2 - x1)[0] * x2[0] + x2[1]);
                        x03 = ((x1 - x2)[2] / (x2 - x1)[0] * x2[0] + x2[2]);
                    }
                    else if (std::abs(dir[1]) > 0.99) // x02 = 0
                    {
                        // std::cout << "parallel to y" << std::endl;
                        x01 = ((x1 - x2)[0] / (x2 - x1)[1] * x2[1] + x2[0]);
                        x03 = ((x1 - x2)[2] / (x2 - x1)[1] * x2[1] + x2[2]);
                    }
                    else if (std::abs(dir[2]) > 0.99) // x03 = 0
                    {
                        // std::cout << "parallel to z" << std::endl;
                        x01 = ((x1 - x2)[0] / (x2 - x1)[2] * x2[2] + x2[0]);
                        x02 = ((x1 - x2)[1] / (x2 - x1)[2] * x2[2] + x2[1]);
                    }
                    else if (std::abs(dir[0]) < 0.1)
                    {
                        // std::cout << "ortho to x" << std::endl;
                        if (exponents[1] < exponents[2]) // x03 = 0
                        {
                            x01 = ((x1 - x2)[0] / (x2 - x1)[2] * x2[2] + x2[0]);
                            x02 = ((x1 - x2)[1] / (x2 - x1)[2] * x2[2] + x2[1]);
                        }
                        else // x02 = 0
                        {
                            x01 = ((x1 - x2)[0] / (x2 - x1)[1] * x2[1] + x2[0]);
                            x03 = ((x1 - x2)[2] / (x2 - x1)[1] * x2[1] + x2[2]);
                        }
                    }
                    else if (std::abs(dir[1]) < 0.1)
                    {
                        // std::cout << "ortho to y" << std::endl;
                        if (exponents[0] < exponents[2]) // x03 = 0
                        {
                            x01 = ((x1 - x2)[0] / (x2 - x1)[2] * x2[2] + x2[0]);
                            x02 = ((x1 - x2)[1] / (x2 - x1)[2] * x2[2] + x2[1]);
                        }
                        else // x01 = 0
                        {
                            x02 = ((x1 - x2)[1] / (x2 - x1)[0] * x2[0] + x2[1]);
                            x03 = ((x1 - x2)[2] / (x2 - x1)[0] * x2[0] + x2[2]);
                        }
                    }
                    else if (std::abs(dir[2]) < 0.1)
                    {
                        // std::cout << "ortho to z" << std::endl;
                        if (exponents[0] < exponents[1]) // x02 = 0
                        {
                            x01 = ((x1 - x2)[0] / (x2 - x1)[1] * x2[1] + x2[0]);
                            x03 = ((x1 - x2)[2] / (x2 - x1)[1] * x2[1] + x2[2]);
                        }
                        else // x01 = 0
                        {
                            x02 = ((x1 - x2)[1] / (x2 - x1)[0] * x2[0] + x2[1]);
                            x03 = ((x1 - x2)[2] / (x2 - x1)[0] * x2[0] + x2[2]);
                        }
                    }
                    else
                    {
                        if (exponents[0] < exponents[1])
                        {
                            if (exponents[1] < exponents[2]) // x03 = 0
                            {
                                // std::cout << "find xy such that z=0 " << (x2 - x1) << std::endl;
                                x01 = ((x1 - x2)[0] / (x2 - x1)[2] * x2[2] + x2[0]);
                                x02 = ((x1 - x2)[1] / (x2 - x1)[2] * x2[2] + x2[1]);
                            }
                            else // x02 = 0
                            {
                                // std::cout << "find xz such that y=0" << std::endl;
                                x01 = ((x1 - x2)[0] / (x2 - x1)[1] * x2[1] + x2[0]);
                                x03 = ((x1 - x2)[2] / (x2 - x1)[1] * x2[1] + x2[2]);
                            }
                        }
                        else
                        {
                            if (exponents[0] < exponents[2]) // x03 = 0
                            {
                                // std::cout << "find xy such that z=0 " << (x2 - x1) << std::endl;
                                x01 = ((x1 - x2)[0] / (x2 - x1)[2] * x2[2] + x2[0]);
                                x02 = ((x1 - x2)[1] / (x2 - x1)[2] * x2[2] + x2[1]);
                            }
                            else // x01 = 0
                            {
                                // std::cout << "find yz such that x=0" << std::endl;
                                x02 = ((x1 - x2)[1] / (x2 - x1)[0] * x2[0] + x2[1]);
                                x03 = ((x1 - x2)[2] / (x2 - x1)[0] * x2[0] + x2[2]);
                            }
                        }
                    }
                    Point3D x0(x01, x02, x03);
                    // std::cout << x0 << std::endl;
                    //  std::cout<<"x0: "<<x0<<std::endl;
                    I += (x0 - x1).dot(dir) * monomial.evaluate(x1);
                    // count++;
                    I += (x2 - x0).dot(dir) * monomial.evaluate(x2);
                    // count++;
                    for (std::size_t i = 0; i < d; i++)
                    {
                        if ((exponents[i] > 0) && (x0[i] != 0.0))
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
                    auto exponents = monomial.getExponents();
                    // std::cout<<"dir: "<<dir<<std::endl;
                    // std::cout<<"abs direction is 1?: "<<(std::abs(dir[0])>0.95)<<" "<<(std::abs(dir[1])>0.95)<<std::endl;
                    if (std::abs(dir[0]) > 0.99)
                        x02 = ((x1 - x2)[1] / (x2 - x1)[0] * x2[0] + x2[1]);
                    else if (std::abs(dir[1]) > 0.99)
                        x01 = ((x1 - x2)[0] / (x2 - x1)[1] * x2[1] + x2[0]);
                    else
                    {
                        if (exponents[0] < exponents[1])
                            x01 = ((x1 - x2)[0] / (x2 - x1)[1] * x2[1] + x2[0]);
                        else
                            x02 = ((x1 - x2)[1] / (x2 - x1)[0] * x2[0] + x2[1]);
                    }
                    Point2D x0(x01, x02);
                    // std::cout<<"x0: "<<x0<<std::endl;
                    I += (x0 - x1).dot(dir) * monomial.evaluate(x1);
                    // count++;
                    I += (x2 - x0).dot(dir) * monomial.evaluate(x2);
                    // count++;
                    for (std::size_t i = 0; i < d; i++)
                    {
                        if ((exponents[i] > 0) && (x0[i] != 0.0))
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
                // auto start = std::chrono::high_resolution_clock::now();

                for (std::size_t e = 0; e < E.numEdges(); e++)
                {
                    bi = ((cross(E[e].getDirection(), E.getOutwardNormal())).dot((E[e][0] - O) / h));
                    // std::cout << "N=" << N << ", d=" << d << ", f=" << f << ", bi=" << bi << std::endl;
                    I += bi * integrateMonomial(N - 1, E[e], monomial, O, h, ex, ey);
                }
                I /= (N + monomial.getOrder());
                I *= std::pow(h, d);

                // auto end = std::chrono::high_resolution_clock::now();
                // elapsed_t += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            }
            // std::cout<<"functions evaluations: "<<count<<std::endl;
            // count=1;
            return I;
        }
        else
            throw std::logic_error("Monomial integration fault.");
        return 0.0;
    }

    Point3D getPolygonCentroid(const Polygon3D &F)
    {
        Point3D X_F(integrateMonomial(2, F, Monomial3D(1, 0, 0, 1.0)),
                    integrateMonomial(2, F, Monomial3D(0, 1, 0, 1.0)),
                    integrateMonomial(2, F, Monomial3D(0, 0, 1, 1.0)));
        return (X_F / F.getArea());
    }

    real getPolyhedronVolume(const Polyhedron<Polygon3D> &P)
    {
        return integrateMonomial(3, P, Monomial3D(0, 0, 0, 1.0));
    }

    Point3D getPolyhedronCentroid(const Polyhedron<Polygon3D> &P)
    {
        Point3D X_P(integrateMonomial(3, P, Monomial3D(1, 0, 0, 1.0)),
                    integrateMonomial(3, P, Monomial3D(0, 1, 0, 1.0)),
                    integrateMonomial(3, P, Monomial3D(0, 0, 1, 1.0)));
        return (X_P / getPolyhedronVolume(P));
    }

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
        auto unit_integrals = MonomialsFaceIntegralsCache::getCacheMonomials(F, poly2D.getOrder());
        real I = 0.0;
        for (const auto &monomialPair : poly2D.getPolynomial())
        {
            auto m = monomialPair.second;
            // std::cout<<m<<std::endl;
            // I += integrateMonomial(2, F, m, X_F, h_F, e_x, e_y);
            // std::cout<<"calling cache monomial "<<m;
            auto order = m.getOrder();
            if (m.getCoefficient() != 0.0)
                I += unit_integrals[order * (order + 1) / 2 + m.getExponents()[1]] * m.getCoefficient();
        }
        return I;
    }

    real integrateMonomial3DRestrictedPolynomial2D(const Point3D &X_P, const real &h_P, const Polygon3D &F, const Monomial3D &m3D, const Polynomial<2> &p2D, const Point3D &X_F = Point3D())
    {
        // auto start = std::chrono::high_resolution_clock::now();

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
            // std::cout<<m<<std::endl;
            // I += integrateMonomial(2, F, m, X_F, h_F, e_x, e_y);
            // std::cout<<"calling cache monomial "<<m;
            auto order = m.getOrder();
            if (m.getCoefficient() != 0.0)
                I += unit_integrals[order * (order + 1) / 2 + m.getExponents()[1]] * m.getCoefficient();
        }

        // auto end = std::chrono::high_resolution_clock::now();
        // elapsed_t += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        return I;
    }
}

namespace Gauss
{
    using traits::real;
    using namespace geometry;

    constexpr int MaxN3 = 216; // N^3

    struct GaussData
    {
        unsigned int N;
        std::array<real, MaxN3> x;
        std::array<real, MaxN3> y;
        std::array<real, MaxN3> z;
        std::array<real, MaxN3> w;
    };

    constexpr GaussData data[] = {
        {1,
         {0.250000000000000},
         {0.250000000000000},
         {0.250000000000000},
         {0.166666666666667}},
        {2,
         {0.544151844011225,
          0.544151844011225,
          0.122514822655441,
          0.122514822655441,
          0.544151844011225,
          0.544151844011225,
          0.122514822655441,
          0.122514822655441},
         {0.293998800631623,
          0.0706797241593969,
          0.565933165072801,
          0.136054976802846,
          0.293998800631623,
          0.0706797241593969,
          0.565933165072801,
          0.136054976802846},
         {0.0342027932367664,
          0.0813956670146703,
          0.0658386870600444,
          0.156682637336818,
          0.127646562120385,
          0.303772764814708,
          0.245713325211713,
          0.584747563204894},
         {0.00916942992147974,
          0.0160270405984766,
          0.0211570064545241,
          0.0369798563588529,
          0.00916942992147974,
          0.0160270405984766,
          0.0211570064545241,
          0.0369798563588529}},
        {3,
         {0.705002209888499,
          0.705002209888499,
          0.705002209888499,
          0.347003766038352,
          0.347003766038352,
          0.347003766038352,
          0.0729940240731496,
          0.0729940240731496,
          0.0729940240731496,
          0.705002209888499,
          0.705002209888499,
          0.705002209888499,
          0.347003766038352,
          0.347003766038352,
          0.347003766038352,
          0.0729940240731496,
          0.0729940240731496,
          0.0729940240731496,
          0.705002209888499,
          0.705002209888499,
          0.705002209888499,
          0.347003766038352,
          0.347003766038352,
          0.347003766038352,
          0.0729940240731496,
          0.0729940240731496,
          0.0729940240731496},
         {0.232357800579865,
          0.120791820133902,
          0.0261332522867348,
          0.514338662174092,
          0.267380320411885,
          0.0578476039361426,
          0.730165028047632,
          0.379578230280591,
          0.0821215678634424,
          0.232357800579865,
          0.120791820133902,
          0.0261332522867348,
          0.514338662174092,
          0.267380320411885,
          0.0578476039361426,
          0.730165028047632,
          0.379578230280591,
          0.0821215678634424,
          0.232357800579865,
          0.120791820133902,
          0.0261332522867348,
          0.514338662174092,
          0.267380320411885,
          0.0578476039361426,
          0.730165028047632,
          0.379578230280591,
          0.0821215678634424},
         {0.00705963113955479,
          0.0196333029354845,
          0.0303014811742758,
          0.0156269392579016,
          0.0434595556538025,
          0.0670742417520586,
          0.0221843026408197,
          0.0616960186091465,
          0.0952198798417150,
          0.0313199947658184,
          0.0871029849887995,
          0.134432268912383,
          0.0693287858937781,
          0.192807956774882,
          0.297574315012753,
          0.0984204739396093,
          0.273713872823130,
          0.422442204031704,
          0.0555803583920821,
          0.154572667042115,
          0.238563056650491,
          0.123030632529655,
          0.342156357895961,
          0.528074388273447,
          0.174656645238399,
          0.485731727037113,
          0.749664528221693},
         {0.000580935315837386,
          0.00190720341498179,
          0.00167168113148370,
          0.00283664869563093,
          0.00931268237947047,
          0.00816265076654669,
          0.00304787709051819,
          0.0100061425721761,
          0.00877047492965105,
          0.000929496505339817,
          0.00305152546397086,
          0.00267468981037392,
          0.00453863791300948,
          0.0149002918071527,
          0.0130602412264747,
          0.00487660334482911,
          0.0160098281154818,
          0.0140327598874417,
          0.000580935315837386,
          0.00190720341498179,
          0.00167168113148370,
          0.00283664869563093,
          0.00931268237947047,
          0.00816265076654669,
          0.00304787709051819,
          0.0100061425721761,
          0.00877047492965105}},
        {4,
         {0.795851417896773,
          0.795851417896773,
          0.795851417896773,
          0.795851417896773,
          0.517047295104368,
          0.517047295104368,
          0.517047295104368,
          0.517047295104368,
          0.238600737551862,
          0.238600737551862,
          0.238600737551862,
          0.238600737551862,
          0.0485005494469972,
          0.0485005494469972,
          0.0485005494469972,
          0.0485005494469972,
          0.795851417896773,
          0.795851417896773,
          0.795851417896773,
          0.795851417896773,
          0.517047295104368,
          0.517047295104368,
          0.517047295104368,
          0.517047295104368,
          0.238600737551862,
          0.238600737551862,
          0.238600737551862,
          0.238600737551862,
          0.0485005494469972,
          0.0485005494469972,
          0.0485005494469972,
          0.0485005494469972,
          0.795851417896773,
          0.795851417896773,
          0.795851417896773,
          0.795851417896773,
          0.517047295104368,
          0.517047295104368,
          0.517047295104368,
          0.517047295104368,
          0.238600737551862,
          0.238600737551862,
          0.238600737551862,
          0.238600737551862,
          0.0485005494469972,
          0.0485005494469972,
          0.0485005494469972,
          0.0485005494469972,
          0.795851417896773,
          0.795851417896773,
          0.795851417896773,
          0.795851417896773,
          0.517047295104368,
          0.517047295104368,
          0.517047295104368,
          0.517047295104368,
          0.238600737551862,
          0.238600737551862,
          0.238600737551862,
          0.238600737551862,
          0.0485005494469972,
          0.0485005494469972,
          0.0485005494469972,
          0.0485005494469972},
         {0.175616803962505,
          0.119139159297124,
          0.0565171086994074,
          0.0116577406689234,
          0.415455300374957,
          0.281846577863780,
          0.133702082267990,
          0.0275786259743970,
          0.654986204816932,
          0.444345324777483,
          0.210788066397987,
          0.0434790928042877,
          0.818518016420534,
          0.555285975747014,
          0.263415975366112,
          0.0543346112272346,
          0.175616803962505,
          0.119139159297124,
          0.0565171086994074,
          0.0116577406689234,
          0.415455300374957,
          0.281846577863780,
          0.133702082267990,
          0.0275786259743970,
          0.654986204816932,
          0.444345324777483,
          0.210788066397987,
          0.0434790928042877,
          0.818518016420534,
          0.555285975747014,
          0.263415975366112,
          0.0543346112272346,
          0.175616803962505,
          0.119139159297124,
          0.0565171086994074,
          0.0116577406689234,
          0.415455300374957,
          0.281846577863780,
          0.133702082267990,
          0.0275786259743970,
          0.654986204816932,
          0.444345324777483,
          0.210788066397987,
          0.0434790928042877,
          0.818518016420534,
          0.555285975747014,
          0.263415975366112,
          0.0543346112272346,
          0.175616803962505,
          0.119139159297124,
          0.0565171086994074,
          0.0116577406689234,
          0.415455300374957,
          0.281846577863780,
          0.133702082267990,
          0.0275786259743970,
          0.654986204816932,
          0.444345324777483,
          0.210788066397987,
          0.0434790928042877,
          0.818518016420534,
          0.555285975747014,
          0.263415975366112,
          0.0543346112272346},
         {0.00198101397470043,
          0.00590236100005810,
          0.0102503254608295,
          0.0133649941129659,
          0.00468646927478463,
          0.0139631692803390,
          0.0242491148180740,
          0.0316174621017319,
          0.00738845483861198,
          0.0220136396042882,
          0.0382299507805671,
          0.0498465213688843,
          0.00923314621657363,
          0.0275098322538483,
          0.0477749046478169,
          0.0622918093484527,
          0.00941575721655393,
          0.0280539152629691,
          0.0487197855050096,
          0.0635238021414710,
          0.0222747832462335,
          0.0663669280461273,
          0.115256015737018,
          0.150277762174051,
          0.0351173176233467,
          0.104630804534349,
          0.181706913503757,
          0.236920460578858,
          0.0438851336893508,
          0.130754202079533,
          0.227074068609678,
          0.296072900492077,
          0.0191160209241683,
          0.0569555075431344,
          0.0989116878988102,
          0.128967039292833,
          0.0452226212744420,
          0.134739198985725,
          0.233994606890624,
          0.305096316747185,
          0.0712957400078597,
          0.212423133136306,
          0.368904282546393,
          0.480999709064992,
          0.0890963004431186,
          0.265459272726456,
          0.461009406577212,
          0.601091938833692,
          0.0265507641660217,
          0.0791070618060454,
          0.137381147942990,
          0.179125847321338,
          0.0628109352458909,
          0.187142957751513,
          0.325001507809568,
          0.423756616819504,
          0.0990246027925943,
          0.295040298066366,
          0.512381245269583,
          0.668073648274966,
          0.123748287915896,
          0.368703642552141,
          0.640308570539074,
          0.834873029977316},
         {5.61425402669510e-05,
          0.000233795515279108,
          0.000366345798555433,
          0.000243985421620605,
          0.000372217075256264,
          0.00155003109035391,
          0.00242882065938498,
          0.00161758872343451,
          0.000778009425931695,
          0.00323988037881461,
          0.00507672939399184,
          0.00338108957856493,
          0.000601372928720176,
          0.00250430944300902,
          0.00392412678076308,
          0.00261345900750740,
          0.000105253918778391,
          0.000438311021534327,
          0.000686811297504770,
          0.000457414673939930,
          0.000697818545806260,
          0.00290593987575818,
          0.00455346144286728,
          0.00303259438036939,
          0.00145858275269461,
          0.00607400564032184,
          0.00951766095289490,
          0.00633873932658917,
          0.00112743130421366,
          0.00469498496963442,
          0.00735680500908297,
          0.00489961445988875,
          0.000105253918778391,
          0.000438311021534327,
          0.000686811297504770,
          0.000457414673939929,
          0.000697818545806260,
          0.00290593987575818,
          0.00455346144286728,
          0.00303259438036939,
          0.00145858275269461,
          0.00607400564032184,
          0.00951766095289490,
          0.00633873932658916,
          0.00112743130421366,
          0.00469498496963442,
          0.00735680500908297,
          0.00489961445988875,
          5.61425402669510e-05,
          0.000233795515279108,
          0.000366345798555433,
          0.000243985421620605,
          0.000372217075256264,
          0.00155003109035392,
          0.00242882065938498,
          0.00161758872343451,
          0.000778009425931696,
          0.00323988037881461,
          0.00507672939399184,
          0.00338108957856493,
          0.000601372928720176,
          0.00250430944300902,
          0.00392412678076308,
          0.00261345900750741}},
        {5,
         {0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150,
          0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150,
          0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150,
          0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150,
          0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.851054212947016,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.634333472630887,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.389886387065519,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.173480320771696,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150,
          0.0345789399182150},
         {0.134269401146344,
          0.103586473561889,
          0.0652345028216781,
          0.0294932643722359,
          0.00592951049099780,
          0.329635544721039,
          0.254308005746508,
          0.160152727938308,
          0.0724068788863313,
          0.0145571321830714,
          0.549996015736950,
          0.424312220482640,
          0.267214393854326,
          0.120810681788372,
          0.0242885357160769,
          0.745078491721125,
          0.574814908126993,
          0.361994799675747,
          0.163661986623795,
          0.0329036302803047,
          0.870293213094632,
          0.671415856030076,
          0.422830105598150,
          0.191166323793956,
          0.0384332743963334,
          0.134269401146344,
          0.103586473561889,
          0.0652345028216781,
          0.0294932643722359,
          0.00592951049099780,
          0.329635544721039,
          0.254308005746508,
          0.160152727938308,
          0.0724068788863313,
          0.0145571321830714,
          0.549996015736950,
          0.424312220482640,
          0.267214393854326,
          0.120810681788372,
          0.0242885357160769,
          0.745078491721125,
          0.574814908126993,
          0.361994799675747,
          0.163661986623795,
          0.0329036302803047,
          0.870293213094632,
          0.671415856030076,
          0.422830105598150,
          0.191166323793956,
          0.0384332743963334,
          0.134269401146344,
          0.103586473561889,
          0.0652345028216781,
          0.0294932643722359,
          0.00592951049099780,
          0.329635544721039,
          0.254308005746508,
          0.160152727938308,
          0.0724068788863313,
          0.0145571321830714,
          0.549996015736950,
          0.424312220482640,
          0.267214393854326,
          0.120810681788372,
          0.0242885357160769,
          0.745078491721125,
          0.574814908126993,
          0.361994799675747,
          0.163661986623795,
          0.0329036302803047,
          0.870293213094632,
          0.671415856030076,
          0.422830105598150,
          0.191166323793956,
          0.0384332743963334,
          0.134269401146344,
          0.103586473561889,
          0.0652345028216781,
          0.0294932643722359,
          0.00592951049099780,
          0.329635544721039,
          0.254308005746508,
          0.160152727938308,
          0.0724068788863313,
          0.0145571321830714,
          0.549996015736950,
          0.424312220482640,
          0.267214393854326,
          0.120810681788372,
          0.0242885357160769,
          0.745078491721125,
          0.574814908126993,
          0.361994799675747,
          0.163661986623795,
          0.0329036302803047,
          0.870293213094632,
          0.671415856030076,
          0.422830105598150,
          0.191166323793956,
          0.0384332743963334,
          0.134269401146344,
          0.103586473561889,
          0.0652345028216781,
          0.0294932643722359,
          0.00592951049099780,
          0.329635544721039,
          0.254308005746508,
          0.160152727938308,
          0.0724068788863313,
          0.0145571321830714,
          0.549996015736950,
          0.424312220482640,
          0.267214393854326,
          0.120810681788372,
          0.0242885357160769,
          0.745078491721125,
          0.574814908126993,
          0.361994799675747,
          0.163661986623795,
          0.0329036302803047,
          0.870293213094632,
          0.671415856030076,
          0.422830105598150,
          0.191166323793956,
          0.0384332743963334},
         {0.000688470393412269,
          0.00212780888992549,
          0.00392690279162669,
          0.00560352704046150,
          0.00670890455016208,
          0.00169021617151184,
          0.00522383682733774,
          0.00964066816216434,
          0.0137568327003139,
          0.0164705687743685,
          0.00282012111543486,
          0.00871595763232123,
          0.0160854287808060,
          0.0229532381913956,
          0.0274810994988124,
          0.00382041237943087,
          0.0118074902013492,
          0.0217908978824723,
          0.0310947054204484,
          0.0372285899889251,
          0.00446245462992894,
          0.0137918067694830,
          0.0254529834709710,
          0.0363203493206217,
          0.0434850684329930,
          0.00338680125632327,
          0.0104673576243388,
          0.0193176633816068,
          0.0275655026012310,
          0.0330032003938848,
          0.00831470213956799,
          0.0256976876550461,
          0.0474254628170509,
          0.0676741639412116,
          0.0810238806942951,
          0.0138730580546826,
          0.0428765224208113,
          0.0791292565731430,
          0.112914159689587,
          0.135188126023001,
          0.0187938037280005,
          0.0580847383280396,
          0.107196244066483,
          0.152964584084757,
          0.183139081291086,
          0.0219522104240708,
          0.0678462123292524,
          0.125211188776624,
          0.178671161296432,
          0.213916656125506,
          0.00733819295331973,
          0.0226796567455475,
          0.0418556421156528,
          0.0597262613403739,
          0.0715081382809929,
          0.0180154913240372,
          0.0556792608113027,
          0.102756899715403,
          0.146629824241391,
          0.175554697593021,
          0.0300587985987655,
          0.0929006962259203,
          0.171449609540077,
          0.244651465573054,
          0.292912538609202,
          0.0407205937535897,
          0.125852385550656,
          0.232262439776279,
          0.331428846302255,
          0.396808024474000,
          0.0475639234935763,
          0.147002602025855,
          0.271295477241817,
          0.387127368143914,
          0.463493892842726,
          0.0112895846503162,
          0.0348919558667561,
          0.0643936208496987,
          0.0918870200795167,
          0.110013076168101,
          0.0277162805085065,
          0.0856608339675593,
          0.158088336613754,
          0.225585484541570,
          0.270085514491747,
          0.0462445391428484,
          0.142924870031029,
          0.263769962507011,
          0.376388771456521,
          0.450636951195403,
          0.0626473837791790,
          0.193620032773272,
          0.357328635486074,
          0.509893108519752,
          0.610476967656913,
          0.0731756365630819,
          0.226158991722457,
          0.417379765707011,
          0.595583574991397,
          0.713071129559946,
          0.0139879155132272,
          0.0432315046011695,
          0.0797843814396788,
          0.113848995640286,
          0.136307372011824,
          0.0343407664765626,
          0.106134684795268,
          0.195873131268641,
          0.279502815782468,
          0.334638826411673,
          0.0572974760820962,
          0.177085434819519,
          0.326813790299348,
          0.466349692954713,
          0.558343977719591,
          0.0776207751277486,
          0.239897280899962,
          0.442733981670085,
          0.631762987184061,
          0.756387458959074,
          0.0906653923572237,
          0.280213397282226,
          0.517137971012664,
          0.737934386967207,
          0.883502717252459},
         {7.67455552179800e-06,
          3.60185932012982e-05,
          7.13399262170556e-05,
          8.14705363128843e-05,
          4.71653365059364e-05,
          5.98013953892923e-05,
          0.000280662785913663,
          0.000555892406098534,
          0.000634831781565256,
          0.000367520038007325,
          0.000166407554052789,
          0.000780991938624512,
          0.00154686516950306,
          0.00176652740822440,
          0.00102268701578053,
          0.000235430746830113,
          0.00110493490770460,
          0.00218848010941898,
          0.00249925473264392,
          0.00144688123847004,
          0.000152536470498619,
          0.000715891501943869,
          0.00141792453255092,
          0.00161927658526933,
          0.000937439821766992,
          1.55037800172007e-05,
          7.27630862707135e-05,
          0.000144117599953650,
          0.000164582987156812,
          9.52812185081395e-05,
          0.000120807996789372,
          0.000566981902660168,
          0.00112298797668545,
          0.00128245763045955,
          0.000742446882427908,
          0.000336168798819303,
          0.00157772357985428,
          0.00312490504969683,
          0.00356865648488400,
          0.00206598473020027,
          0.000475606241660781,
          0.00223213809499741,
          0.00442106570107948,
          0.00504887813656487,
          0.00292292216383615,
          0.000308147081155882,
          0.00144621070637858,
          0.00286442517370849,
          0.00327118722298826,
          0.00189377231486030,
          1.84274965775891e-05,
          8.64848134932766e-05,
          0.000171295424533232,
          0.000195620192572181,
          0.000113249435042247,
          0.000143590075769373,
          0.000673903851785407,
          0.00133476204345559,
          0.00152430462570917,
          0.000882458172768385,
          0.000399563808494583,
          0.00187525208922537,
          0.00371420241029557,
          0.00424163688396196,
          0.00245558995953754,
          0.000565296487744314,
          0.00265307667295563,
          0.00525479418474412,
          0.00600100004508526,
          0.00347412941301362,
          0.000366257730507927,
          0.00171893840164767,
          0.00340460100870314,
          0.00388807060532281,
          0.00225090157446145,
          1.55037800172007e-05,
          7.27630862707136e-05,
          0.000144117599953650,
          0.000164582987156812,
          9.52812185081396e-05,
          0.000120807996789372,
          0.000566981902660169,
          0.00112298797668545,
          0.00128245763045955,
          0.000742446882427908,
          0.000336168798819303,
          0.00157772357985428,
          0.00312490504969683,
          0.00356865648488400,
          0.00206598473020028,
          0.000475606241660782,
          0.00223213809499741,
          0.00442106570107948,
          0.00504887813656487,
          0.00292292216383615,
          0.000308147081155882,
          0.00144621070637859,
          0.00286442517370849,
          0.00327118722298826,
          0.00189377231486030,
          7.67455552179800e-06,
          3.60185932012982e-05,
          7.13399262170555e-05,
          8.14705363128843e-05,
          4.71653365059364e-05,
          5.98013953892923e-05,
          0.000280662785913663,
          0.000555892406098534,
          0.000634831781565255,
          0.000367520038007325,
          0.000166407554052789,
          0.000780991938624512,
          0.00154686516950306,
          0.00176652740822439,
          0.00102268701578053,
          0.000235430746830113,
          0.00110493490770460,
          0.00218848010941898,
          0.00249925473264392,
          0.00144688123847004,
          0.000152536470498619,
          0.000715891501943869,
          0.00141792453255092,
          0.00161927658526933,
          0.000937439821766992}},
        {6,
         {0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.886805616177562,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.715681127311714,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.509036413164752,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.302436918022891,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.131563941657985,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667,
          0.025904555093667},
         {0.104925044101914,
          0.087072953027846,
          0.063238478326491,
          0.038144774373970,
          0.016761665846764,
          0.003318454908227,
          0.263548148312915,
          0.218707704486042,
          0.158840856420180,
          0.095811107253949,
          0.042101540527004,
          0.008335213521318,
          0.455096571592547,
          0.377665816017406,
          0.274287372718776,
          0.165447212249001,
          0.072701200426780,
          0.014393298231562,
          0.646603079311139,
          0.536589143558839,
          0.389708626447115,
          0.235068079131623,
          0.103294164359602,
          0.020450057282137,
          0.804993045098109,
          0.668030423084586,
          0.485170491669912,
          0.292649656149370,
          0.128596795420904,
          0.025459442447314,
          0.902933556108203,
          0.749307201071183,
          0.544199381634270,
          0.328255251806044,
          0.144242689425179,
          0.028556998157275,
          0.104925044101914,
          0.087072953027846,
          0.063238478326491,
          0.038144774373970,
          0.016761665846764,
          0.003318454908227,
          0.263548148312915,
          0.218707704486042,
          0.158840856420180,
          0.095811107253949,
          0.042101540527004,
          0.008335213521318,
          0.455096571592547,
          0.377665816017406,
          0.274287372718776,
          0.165447212249001,
          0.072701200426780,
          0.014393298231562,
          0.646603079311139,
          0.536589143558839,
          0.389708626447115,
          0.235068079131623,
          0.103294164359602,
          0.020450057282137,
          0.804993045098109,
          0.668030423084586,
          0.485170491669912,
          0.292649656149370,
          0.128596795420904,
          0.025459442447314,
          0.902933556108203,
          0.749307201071183,
          0.544199381634270,
          0.328255251806044,
          0.144242689425179,
          0.028556998157275,
          0.104925044101914,
          0.087072953027846,
          0.063238478326491,
          0.038144774373970,
          0.016761665846764,
          0.003318454908227,
          0.263548148312915,
          0.218707704486042,
          0.158840856420180,
          0.095811107253949,
          0.042101540527004,
          0.008335213521318,
          0.455096571592547,
          0.377665816017406,
          0.274287372718776,
          0.165447212249001,
          0.072701200426780,
          0.014393298231562,
          0.646603079311139,
          0.536589143558839,
          0.389708626447115,
          0.235068079131623,
          0.103294164359602,
          0.020450057282137,
          0.804993045098109,
          0.668030423084586,
          0.485170491669912,
          0.292649656149370,
          0.128596795420904,
          0.025459442447314,
          0.902933556108203,
          0.749307201071183,
          0.544199381634270,
          0.328255251806044,
          0.144242689425179,
          0.028556998157275,
          0.104925044101914,
          0.087072953027846,
          0.063238478326491,
          0.038144774373970,
          0.016761665846764,
          0.003318454908227,
          0.263548148312915,
          0.218707704486042,
          0.158840856420180,
          0.095811107253949,
          0.042101540527004,
          0.008335213521318,
          0.455096571592547,
          0.377665816017406,
          0.274287372718776,
          0.165447212249001,
          0.072701200426780,
          0.014393298231562,
          0.646603079311139,
          0.536589143558839,
          0.389708626447115,
          0.235068079131623,
          0.103294164359602,
          0.020450057282137,
          0.804993045098109,
          0.668030423084586,
          0.485170491669912,
          0.292649656149370,
          0.128596795420904,
          0.025459442447314,
          0.902933556108203,
          0.749307201071183,
          0.544199381634270,
          0.328255251806044,
          0.144242689425179,
          0.028556998157275,
          0.104925044101914,
          0.087072953027846,
          0.063238478326491,
          0.038144774373970,
          0.016761665846764,
          0.003318454908227,
          0.263548148312915,
          0.218707704486042,
          0.158840856420180,
          0.095811107253949,
          0.042101540527004,
          0.008335213521318,
          0.455096571592547,
          0.377665816017406,
          0.274287372718776,
          0.165447212249001,
          0.072701200426780,
          0.014393298231562,
          0.646603079311139,
          0.536589143558839,
          0.389708626447115,
          0.235068079131623,
          0.103294164359602,
          0.020450057282137,
          0.804993045098109,
          0.668030423084586,
          0.485170491669912,
          0.292649656149370,
          0.128596795420904,
          0.025459442447314,
          0.902933556108203,
          0.749307201071183,
          0.544199381634270,
          0.328255251806044,
          0.144242689425179,
          0.028556998157275,
          0.104925044101914,
          0.087072953027846,
          0.063238478326491,
          0.038144774373970,
          0.016761665846764,
          0.003318454908227,
          0.263548148312915,
          0.218707704486042,
          0.158840856420180,
          0.095811107253949,
          0.042101540527004,
          0.008335213521318,
          0.455096571592547,
          0.377665816017406,
          0.274287372718776,
          0.165447212249001,
          0.072701200426780,
          0.014393298231562,
          0.646603079311139,
          0.536589143558839,
          0.389708626447115,
          0.235068079131623,
          0.103294164359602,
          0.020450057282137,
          0.804993045098109,
          0.668030423084586,
          0.485170491669912,
          0.292649656149370,
          0.128596795420904,
          0.025459442447314,
          0.902933556108203,
          0.749307201071183,
          0.544199381634270,
          0.328255251806044,
          0.144242689425179,
          0.028556998157275},
         {0.000279216264273,
          0.000881996455634,
          0.001686773283281,
          0.002534068292459,
          0.003256074145804,
          0.003709987428478,
          0.000701328553711,
          0.002215377031198,
          0.004236795697705,
          0.006365010488130,
          0.008178527054634,
          0.009318655287769,
          0.001211058481711,
          0.003825526751514,
          0.007316124999954,
          0.010991139455319,
          0.014122731072357,
          0.016091511552876,
          0.001720676868121,
          0.005435324131009,
          0.010394780468332,
          0.015616256110385,
          0.020065634350339,
          0.022862885748512,
          0.002142168752388,
          0.006766744952680,
          0.012941054953286,
          0.019441567727644,
          0.024980852418322,
          0.028463310193376,
          0.002402798460380,
          0.007590029653804,
          0.014515544997458,
          0.021806950993596,
          0.028020179858748,
          0.031926335324281,
          0.001400787338718,
          0.004424847782639,
          0.008462295936303,
          0.012713051615257,
          0.016335249843852,
          0.018612466684717,
          0.003518463227336,
          0.011114223964952,
          0.021255387058234,
          0.031932330753686,
          0.041030479285713,
          0.046750336607231,
          0.006075704049849,
          0.019192110643691,
          0.036703933759343,
          0.055140946130676,
          0.070851685254705,
          0.080728770233996,
          0.008632385284442,
          0.027268229679834,
          0.052149099934050,
          0.078344482885154,
          0.100666365601832,
          0.114699768534046,
          0.010746948690666,
          0.033947774062241,
          0.064923388239598,
          0.097535514231614,
          0.125325298900695,
          0.142796282446779,
          0.012054489983069,
          0.038078073522041,
          0.072822375515840,
          0.109402297632597,
          0.140573162086061,
          0.160169775246924,
          0.003148058303483,
          0.009944178119529,
          0.019017733993228,
          0.028570666363006,
          0.036711010650264,
          0.041828712093283,
          0.007907215515281,
          0.024977542323929,
          0.047768277077438,
          0.071763097938017,
          0.092209814752857,
          0.105064331522142,
          0.013654228629127,
          0.043131374480124,
          0.082486556130206,
          0.123920961112862,
          0.159228478097232,
          0.181425737112820,
          0.019399984153509,
          0.061281234126148,
          0.117197237959671,
          0.176067410849475,
          0.226232476090537,
          0.257770432927963,
          0.024152146530490,
          0.076292502842908,
          0.145905524549545,
          0.219196359771832,
          0.281649710085242,
          0.320913110961370,
          0.027090648406489,
          0.085574728025068,
          0.163657307276856,
          0.245865165941442,
          0.315916983199671,
          0.359957416087714,
          0.005121281417040,
          0.016177252675063,
          0.030938171502719,
          0.046478943085461,
          0.059721707325410,
          0.068047216820928,
          0.012863508860090,
          0.040633625878316,
          0.077709739190668,
          0.116744667496320,
          0.150007517408426,
          0.170919327644826,
          0.022212786613574,
          0.070166396337718,
          0.134189657986267,
          0.201595413473385,
          0.259033908311236,
          0.295144551490866,
          0.031560018512461,
          0.099692704292122,
          0.190657217570323,
          0.286427591996010,
          0.368036441526969,
          0.419342591767009,
          0.039290866713416,
          0.124113132414521,
          0.237360042122558,
          0.356590042420813,
          0.458189552835869,
          0.522063504933331,
          0.044071240391642,
          0.139213515810082,
          0.266238755995207,
          0.399975027158847,
          0.513935772281482,
          0.585581030661344,
          0.006868552381806,
          0.021696583011952,
          0.041493609559644,
          0.062336557833211,
          0.080097468131822,
          0.091263462229494,
          0.017252261148035,
          0.054496944237293,
          0.104222629209872,
          0.156575434680651,
          0.201186852875570,
          0.229233322559737,
          0.029791311192852,
          0.094105660174152,
          0.179972280357129,
          0.270375428455571,
          0.347410701153764,
          0.395841518369690,
          0.042327617381528,
          0.133705708738436,
          0.255705355595944,
          0.384150519960331,
          0.493602552015674,
          0.562413256160925,
          0.052696064553240,
          0.166457861195188,
          0.318342178432505,
          0.478250887961030,
          0.614513964020416,
          0.700180333447922,
          0.059107398815061,
          0.186710170313110,
          0.357073687756223,
          0.536437895467693,
          0.689279593395093,
          0.785368671502134,
          0.007990123456251,
          0.025239434338958,
          0.048269132212666,
          0.072515541156008,
          0.093176643829871,
          0.106165941485733,
          0.020069395821660,
          0.063395791171046,
          0.121241220570401,
          0.182142754946208,
          0.234038805106649,
          0.266665003879199,
          0.034655956760990,
          0.109472244066328,
          0.209360089116519,
          0.314525235130928,
          0.404139655336112,
          0.460478777050810,
          0.049239325797849,
          0.155538614287261,
          0.297459675061662,
          0.446878746735100,
          0.574203283267168,
          0.654250138946459,
          0.061300844491518,
          0.193638890304749,
          0.370324511718817,
          0.556344834465000,
          0.714858410502789,
          0.814513305701324,
          0.068759090337751,
          0.217198214181346,
          0.415380518274606,
          0.624033242106693,
          0.801832575622406,
          0.913612111424778},
         {0.000001370638174,
          0.000006894546463,
          0.000015475402729,
          0.000022083884603,
          0.000021260391882,
          0.000011342169429,
          0.000011767311518,
          0.000059191606893,
          0.000132860654397,
          0.000189596316906,
          0.000182526401914,
          0.000097375692186,
          0.000038392426938,
          0.000193120530504,
          0.000433475646412,
          0.000618583330038,
          0.000595516787239,
          0.000317701213401,
          0.000070795049603,
          0.000356111312226,
          0.000799322479379,
          0.001140658224197,
          0.001098123871152,
          0.000585836191028,
          0.000080375880130,
          0.000404304542544,
          0.000907496331276,
          0.001295025700393,
          0.001246735091362,
          0.000665118532009,
          0.000046812880715,
          0.000235476865600,
          0.000528548084788,
          0.000754254678619,
          0.000726129045449,
          0.000387381319499,
          0.000002886181520,
          0.000014517990936,
          0.000032586879753,
          0.000046502498491,
          0.000044768452617,
          0.000023883443797,
          0.000024778674403,
          0.000124641006768,
          0.000279767463548,
          0.000399236936769,
          0.000384349668649,
          0.000205046035169,
          0.000080843737769,
          0.000406657947180,
          0.000912778750463,
          0.001302563878094,
          0.001253992175651,
          0.000668990101264,
          0.000149074619187,
          0.000749871051096,
          0.001683150091555,
          0.002401907920794,
          0.002312342442444,
          0.001233607541391,
          0.000169249174756,
          0.000851352545880,
          0.001910934037860,
          0.002726962749597,
          0.002625276202417,
          0.001400554027851,
          0.000098574863706,
          0.000495848569406,
          0.001112974775829,
          0.001588249879271,
          0.001529025144241,
          0.000815716960555,
          0.000003743426828,
          0.000018830082715,
          0.000042265740755,
          0.000060314536420,
          0.000058065449243,
          0.000030977235363,
          0.000032138364780,
          0.000161661518972,
          0.000362863188360,
          0.000517817139814,
          0.000498508098271,
          0.000265948216918,
          0.000104855711503,
          0.000527442317325,
          0.001183889661293,
          0.001689447642889,
          0.001626449313531,
          0.000867691610911,
          0.000193352332455,
          0.000972595587096,
          0.002183074475371,
          0.003115315681224,
          0.002999147722918,
          0.001600010094018,
          0.000219519076309,
          0.001104218822658,
          0.002478514152184,
          0.003536917357400,
          0.003405028165372,
          0.001816542544197,
          0.000127853285307,
          0.000643124080824,
          0.001443547332489,
          0.002059987275849,
          0.001983171780900,
          0.001057998858598,
          0.000003743426828,
          0.000018830082715,
          0.000042265740755,
          0.000060314536420,
          0.000058065449243,
          0.000030977235363,
          0.000032138364780,
          0.000161661518972,
          0.000362863188360,
          0.000517817139814,
          0.000498508098271,
          0.000265948216918,
          0.000104855711503,
          0.000527442317325,
          0.001183889661293,
          0.001689447642889,
          0.001626449313531,
          0.000867691610911,
          0.000193352332455,
          0.000972595587096,
          0.002183074475371,
          0.003115315681224,
          0.002999147722918,
          0.001600010094018,
          0.000219519076309,
          0.001104218822658,
          0.002478514152184,
          0.003536917357400,
          0.003405028165372,
          0.001816542544197,
          0.000127853285307,
          0.000643124080824,
          0.001443547332489,
          0.002059987275849,
          0.001983171780900,
          0.001057998858598,
          0.000002886181520,
          0.000014517990936,
          0.000032586879753,
          0.000046502498491,
          0.000044768452617,
          0.000023883443797,
          0.000024778674403,
          0.000124641006768,
          0.000279767463548,
          0.000399236936769,
          0.000384349668649,
          0.000205046035169,
          0.000080843737769,
          0.000406657947180,
          0.000912778750463,
          0.001302563878094,
          0.001253992175651,
          0.000668990101264,
          0.000149074619187,
          0.000749871051096,
          0.001683150091555,
          0.002401907920794,
          0.002312342442444,
          0.001233607541391,
          0.000169249174756,
          0.000851352545880,
          0.001910934037860,
          0.002726962749597,
          0.002625276202417,
          0.001400554027851,
          0.000098574863706,
          0.000495848569406,
          0.001112974775829,
          0.001588249879271,
          0.001529025144241,
          0.000815716960555,
          0.000001370638174,
          0.000006894546463,
          0.000015475402729,
          0.000022083884603,
          0.000021260391882,
          0.000011342169429,
          0.000011767311518,
          0.000059191606893,
          0.000132860654397,
          0.000189596316906,
          0.000182526401914,
          0.000097375692186,
          0.000038392426938,
          0.000193120530504,
          0.000433475646412,
          0.000618583330038,
          0.000595516787239,
          0.000317701213401,
          0.000070795049603,
          0.000356111312226,
          0.000799322479379,
          0.001140658224197,
          0.001098123871152,
          0.000585836191028,
          0.000080375880130,
          0.000404304542544,
          0.000907496331276,
          0.001295025700393,
          0.001246735091362,
          0.000665118532009,
          0.000046812880715,
          0.000235476865600,
          0.000528548084788,
          0.000754254678619,
          0.000726129045449,
          0.000387381319499}}};

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
        // std::cout << XYZ << std::endl;
        // std::cout << W[0] << std::endl;

        /*
        // Extract the specified column from the matrix
        Eigen::VectorXd Xv = XYZ.col(0);
        Eigen::VectorXd Yv = XYZ.col(1);
        Eigen::VectorXd Zv = XYZ.col(2);

        // Convert Eigen column vector to std::vector
        std::vector<real> X(Xv.data(), Xv.data() + Xv.size());
        std::vector<real> Y(Yv.data(), Yv.data() + Yv.size());
        std::vector<real> Z(Zv.data(), Zv.data() + Zv.size());
        */

        std::vector<std::array<real, 4>> result;
        result.reserve(NP);
        for (Eigen::Index i = 0; i < XYZ.rows(); i++)
        {
            std::array<real, 4> arr = {XYZ(i, 0),
                                       XYZ(i, 1),
                                       XYZ(i, 2),
                                       W[i]};
            result.push_back(arr);
            /*
            for (std::size_t j = 0; j < 4; j++)
            {
                std::cout << arr[j] << " ";
            }
            std::cout << std::endl;
            */
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
    /*
        real integrateFunctionOverTetrahedron(const std::function<real(Point3D)> &func,
                                              const Point3D &v1,
                                              const Point3D &v2,
                                              const Point3D &v3,
                                              const Point3D &v4,
                                              unsigned int N_)
        {
            real I = 0.0;
            auto xyzw = gaussPointsTetrahedron(v1, v2, v3, v4, N_);
            for (std::size_t i = 0; i < xyzw[0].size(); ++i)
            {
                real x = xyzw[0][i];
                real y = xyzw[1][i];
                real z = xyzw[2][i];
                real w = xyzw[3][i];
                I += w * func(Point3D(x, y, z));
            }
            return I;
        }
    */
    /*
     real integrateMonomialFunctionOverPolyhedron(const Monomial3D &monomial,
                                                  const std::function<real(real, real, real)> &func,
                                                  const Polyhedron<Polygon3D> &P,
                                                  const unsigned int &N_)
     {
         real I = 0.0;
         auto xyzw = gaussPointsPolyhedron(P, N_);
         for (std::size_t i = 0; i < xyzw[0].size(); ++i)
         {
             real x = xyzw[0][i];
             real y = xyzw[1][i];
             real z = xyzw[2][i];
             real w = xyzw[3][i];
             I += w * (func(x, y, z) * monomial.evaluate(Point3D(x, y, z)));
         }
         return I;
     }
    */
    real integrateFunctionOverPolyhedron(const std::function<real(real, real, real)> &func,
                                         const Polyhedron<Polygon3D> &P,
                                         const unsigned int &N_)
    {
        real I = 0.0;
        auto xyzw = gaussPointsPolyhedron(P, N_);
        //std::cout<<std::setprecision(16);
        for (std::size_t i = 0; i < xyzw.size(); ++i)
        {
            real x = xyzw[i][0];
            real y = xyzw[i][1];
            real z = xyzw[i][2];
            real w = xyzw[i][3];
            I += w * func(x, y, z);
            //if (P.getId()==0)
              //  std::cout<<xyzw[i][0]<<" "<<xyzw[i][1]<<" "<<xyzw[i][2]<<" "<<xyzw[i][3]<<" "<<std::endl;
        }
        //std::cout<<std::endl;
        return I;
    }
}

#endif // __INTEGRATION_HPP_