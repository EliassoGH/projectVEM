#ifndef __MONOMIAL_HPP_
#define __MONOMIAL_HPP_

#include "traits.hpp"
#include "point.hpp"
#include <vector>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <utility>
#include <cmath>
#include <map>

using namespace traits;
using namespace geometry;

template <unsigned int Dimension>
class Monomial
{
private:
    std::vector<unsigned int> exponents;
    real coefficient;

public:
    Monomial() : exponents(Dimension, 0), coefficient(0.0) {}

    Monomial(std::vector<unsigned int> exponents, real coefficient) : exponents(std::move(exponents)), coefficient(coefficient) {}

    // Getter for the exponents
    const std::vector<unsigned int> &getExponents() const
    {
        return exponents;
    }

    // Getter for the coefficient
    real getCoefficient() const
    {
        return coefficient;
    }

    // Setter for the coefficient
    void setCoefficient(real coeff)
    {
        coefficient = coeff;
    }

    // Compute the product of two monomials
    Monomial<Dimension> operator*(const Monomial<Dimension> &other) const
    {
        std::vector<unsigned int> newExponents(exponents.size(), 0);

        std::transform(exponents.begin(), exponents.end(), other.exponents.begin(), newExponents.begin(), std::plus<unsigned int>());

        return Monomial<Dimension>(std::move(newExponents), coefficient * other.coefficient);
    }

    // Compute the derivative with respect to a variable
    Monomial<Dimension> derivative(unsigned int variableIndex) const
    {
        if (variableIndex < exponents.size())
        {
            if (exponents[variableIndex] > 0)
            {
                std::vector<unsigned int> newExponents = exponents;
                newExponents[variableIndex]--;

                real newCoefficient = coefficient * exponents[variableIndex];

                return Monomial<Dimension>(std::move(newExponents), newCoefficient);
            }
            else
            {
                // The derivative with respect to the variable is zero
                return Monomial<Dimension>(std::vector<unsigned int>(exponents.size(), 0), 0);
            }
        }
        else
        {
            // Invalid variable index
            throw std::out_of_range("Invalid variable index for derivative.");
        }
    }

    // Evaluate the monomial at a point
    template <typename PointType>
    real evaluate(const PointType &point) const
    {
        if (point.getDimension() != exponents.size())
        {
            throw std::invalid_argument("Number of coordinates in the point does not match the dimension of the monomial.");
        }

        real result = coefficient;

        for (unsigned int i = 0; i < exponents.size(); ++i)
        {
            result *= std::pow(point[i], exponents[i]);
        }

        return result;
    }

    // Get the monomial order
    unsigned int getOrder() const
    {
        return std::accumulate(exponents.begin(), exponents.end(), 0);
    }

    friend std::ostream &operator<<(std::ostream &os, const Monomial &monomial)
    {
        // Print the coefficient
        os << monomial.coefficient << "(";

        // Print each variable with its exponent
        for (unsigned int i = 0; i < monomial.exponents.size(); ++i)
        {
            if (monomial.exponents[i] != 0)
            {
                os << "x" << (i + 1);
                if (monomial.exponents[i] != 1)
                {
                    os << "^" << monomial.exponents[i];
                }
            }
        }
        os << ")";
        return os;
    }
};

class Monomial2D : public Monomial<2>
{
private:
    static unsigned int order;
    static std::vector<Monomial2D> monomials_ordered;
    static std::vector<std::pair<std::pair<real, std::size_t>,
                                 std::pair<real, std::size_t>>>
        laplacians_to_monomials_ordered;

public:
    using Monomial<2>::Monomial;

    // Default constructor
    Monomial2D() : Monomial<2>() {}

    // Constructor with exponents and coefficient
    Monomial2D(unsigned int expX, unsigned int expY, double coeff) : Monomial<2>({expX, expY}, coeff) {}

    // Constructor for implicit conversion from Monomial<2> to Monomial2D
    Monomial2D(const Monomial<2> &monomial)
        : Monomial<2>(monomial.getExponents(), monomial.getCoefficient()) {}

    // Fill monomials_ordered
    static void computeMonomialsUpToOrder(unsigned int order_)
    {
        if (order_ > order)
        {
            std::size_t j;
            for (std::size_t kk = (order + 1); kk <= order_; kk++)
            {
                for (std::size_t i = kk; i > 0; i--)
                {
                    j = kk - i;
                    monomials_ordered.emplace_back(i, j, 1.0);
                }
                monomials_ordered.emplace_back(0, kk, 1.0);
            }
            order = order_;
        }
    }

    static const std::vector<Monomial2D> getMonomialsOrdered(unsigned int order_)
    {
        if (order_ > order) // Check if the cache for n doesn't exist
        {
            computeMonomialsUpToOrder(order_); // Initialize the cache for n
        }
        std::size_t endIndex = (order_ + 1) * (order_ + 2) / 2;
        return std::vector<Monomial2D>(monomials_ordered.begin(), monomials_ordered.begin() + endIndex);
    }

    static void computeLaplaciansToMonomialsOrdered(unsigned int order_)
    {
        auto monomials = getMonomialsOrdered(order_);
        std::size_t n_kminus2 = 0;
        if (order_ > 1)
            n_kminus2 = getMonomialsOrdered(order_ - 2).size();
        for (std::size_t i = laplacians_to_monomials_ordered.size(); i < monomials.size(); i++)
        {
            auto mdxx = monomials[i].dx().dx();
            auto mdyy = monomials[i].dy().dy();
            real c1(0.0), c2(0.0);
            std::size_t m1(0), m2(0);
            if (mdxx.getCoefficient() != 0.0)
            {
                for (std::size_t j = 0; j < n_kminus2; j++)
                {
                    if (monomials[j].getExponents() == mdxx.getExponents())
                    {
                        c1 = mdxx.getCoefficient();
                        m1 = j;
                        break;
                    }
                }
            }
            if (mdyy.getCoefficient() != 0.0)
            {
                for (std::size_t j = 0; j < n_kminus2; j++)
                {
                    if (monomials[j].getExponents() == mdyy.getExponents())
                    {
                        c2 = mdyy.getCoefficient();
                        m2 = j;
                        break;
                    }
                }
            }
            laplacians_to_monomials_ordered.emplace_back(std::make_pair(std::make_pair(c1, m1),
                                                                        std::make_pair(c2, m2)));
        }
    }

    static const std::vector<std::pair<std::pair<real, std::size_t>,
                                       std::pair<real, std::size_t>>>
    getLaplaciansToMonomialsOrdered(unsigned int order_)
    {
        if ((laplacians_to_monomials_ordered.size() < ((order_ + 1) * (order_ + 2) / 2))) // Check if the cache for n doesn't exist
        {
            computeLaplaciansToMonomialsOrdered(order_); // Initialize the cache for n
        }
        std::size_t endIndex = (order_ + 1) * (order_ + 2) / 2;
        return std::vector<std::pair<std::pair<real, std::size_t>,
                                     std::pair<real, std::size_t>>>(laplacians_to_monomials_ordered.begin(), laplacians_to_monomials_ordered.begin() + endIndex);
    }

    // Method to compute the product of two monomials and return a new Monomial2D instance
    Monomial2D operator*(const Monomial2D &other) const
    {
        return Monomial::operator*(other);
    }

    // Method to compute the derivative with respect to a variable and return a new Monomial2D instance
    Monomial2D derivative(unsigned int variableIndex) const
    {
        return Monomial::derivative(variableIndex);
    }

    // Compute the derivative with respect to x
    Monomial2D dx() const
    {
        return this->derivative(0);
    }

    // Compute the derivative with respect to y
    Monomial2D dy() const
    {
        return this->derivative(1);
    }

    // Overriding the output stream operator to print "x, y" for Monomial2D
    friend std::ostream &operator<<(std::ostream &os, const Monomial2D &monomial)
    {
        // Print the coefficient
        os << monomial.getCoefficient() << "(";

        // Print variables with their exponents
        if (monomial.getExponents()[0] > 0)
        {
            os << "x";
            if (monomial.getExponents()[0] != 1)
            {
                os << "^" << monomial.getExponents()[0];
            }
        }

        if (monomial.getExponents()[1] > 0)
        {
            if (monomial.getExponents()[0] > 0)
            {
                os << " ";
            }
            os << "y";
            if (monomial.getExponents()[1] != 1)
            {
                os << "^" << monomial.getExponents()[1];
            }
        }

        os << ")";
        return os;
    }
};

unsigned int Monomial2D::order = 0;
std::vector<Monomial2D> Monomial2D::monomials_ordered = {Monomial2D(0, 0, 1.0)};
std::vector<std::pair<std::pair<real, std::size_t>,
                      std::pair<real, std::size_t>>>
    Monomial2D::laplacians_to_monomials_ordered = {};

class Monomial3D : public Monomial<3>
{
private:
    static unsigned int order;
    static std::vector<Monomial3D> monomials_ordered;

public:
    // Default constructor
    Monomial3D() : Monomial<3>() {}

    // Constructor with exponents and coefficient
    Monomial3D(unsigned int expX, unsigned int expY, unsigned int expZ, double coeff) : Monomial<3>({expX, expY, expZ}, coeff) {}

    // Constructor for implicit conversion from Monomial<3> to Monomial3D
    Monomial3D(const Monomial<3> &monomial)
        : Monomial<3>(monomial.getExponents(), monomial.getCoefficient()) {}

    // Fill monomials_ordered
    static void computeMonomialsUpToOrder(unsigned int order_)
    {
        if (order_ > order)
        {
            std::size_t j;
            for (std::size_t kk = (order + 1); kk <= order_; kk++)
            {
                for (std::size_t k = 0; k < (kk / 3 + 1); k++)
                {
                    for (std::size_t i = (kk - 2 * k); i > k; i--)
                    {
                        j = kk - i - k;
                        monomials_ordered.emplace_back(i, j, k, 1.0);
                    }
                    for (std::size_t i = (kk - 2 * k); i > k; i--)
                    {
                        j = kk - i - k;
                        monomials_ordered.emplace_back(k, i, j, 1.0);
                    }
                    for (std::size_t i = (kk - 2 * k); i > k; i--)
                    {
                        j = kk - i - k;
                        monomials_ordered.emplace_back(j, k, i, 1.0);
                    }
                }
                if (kk % 3 == 0)
                {
                    monomials_ordered.emplace_back(kk / 3, kk / 3, kk / 3, 1.0);
                }
            }
        }
        order = order_;
    }

    static const std::vector<Monomial3D>
    getMonomialsOrdered(unsigned int order_)
    {
        if (order_ > order) // Check if the cache for n doesn't exist
        {
            computeMonomialsUpToOrder(order_); // Initialize the cache for n
        }
        std::size_t endIndex = (order_ + 1) * (order_ + 2) * (order_ + 3) / 6;
        return std::vector<Monomial3D>(monomials_ordered.begin(), monomials_ordered.begin() + endIndex);
    }

    // Method to compute the product of two monomials and return a new Monomial2D instance
    Monomial3D operator*(const Monomial3D &other) const
    {
        return Monomial::operator*(other);
    }

    // Method to compute the derivative with respect to a variable and return a new Monomial2D instance
    Monomial3D derivative(unsigned int variableIndex) const
    {
        return Monomial::derivative(variableIndex);
    }

    // Compute the derivative with respect to x
    Monomial3D dx() const
    {
        return this->derivative(0);
    }

    // Compute the derivative with respect to y
    Monomial3D dy() const
    {
        return this->derivative(1);
    }

    // Compute the derivative with respect to z
    Monomial3D dz() const
    {
        return this->derivative(2);
    }

    // Overriding the output stream operator to print "x, y, z" for Monomial3D
    friend std::ostream &operator<<(std::ostream &os, const Monomial3D &monomial)
    {
        // Print the coefficient
        os << monomial.getCoefficient() << "(";

        // Print variables with their exponents
        if (monomial.getExponents()[0] > 0)
        {
            os << "x";
            if (monomial.getExponents()[0] != 1)
            {
                os << "^" << monomial.getExponents()[0];
            }
        }

        if (monomial.getExponents()[1] > 0)
        {
            if (monomial.getExponents()[0] > 0)
            {
                os << " ";
            }
            os << "y";
            if (monomial.getExponents()[1] != 1)
            {
                os << "^" << monomial.getExponents()[1];
            }
        }

        if (monomial.getExponents()[2] > 0)
        {
            if (monomial.getExponents()[0] > 0 || monomial.getExponents()[1] > 0)
            {
                os << " ";
            }
            os << "z";
            if (monomial.getExponents()[2] != 1)
            {
                os << "^" << monomial.getExponents()[2];
            }
        }

        os << ")";
        return os;
    }
};

unsigned int Monomial3D::order = 0;
std::vector<Monomial3D> Monomial3D::monomials_ordered = {Monomial3D(0, 0, 0, 1.0)};

unsigned int factorial(unsigned int n)
{
    return std::tgamma(n + 1); // Factorial using tgamma
}

unsigned int multinomialCoeff(unsigned int n, const std::vector<unsigned int> &k_values)
{
    unsigned int denominator = 1;
    for (unsigned int k : k_values)
    {
        denominator *= factorial(k);
    }

    return factorial(n) / denominator;
}

template <unsigned int Dimension>
class Polynomial
{
protected:
    // Map of monomials
    std::map<std::vector<unsigned int>, Monomial<Dimension>> polynomial;

public:
    Polynomial() = default;

    // Constructor
    Polynomial(const std::vector<Monomial<Dimension>> &monomials)
    {
        for (const Monomial<Dimension> &monomial : monomials)
        {
            addMonomial(monomial);
        }
    }

    // Method to add a monomial to the polynomial
    void addMonomial(const Monomial<Dimension> &monomial)
    {
        auto it = polynomial.find(monomial.getExponents());
        if (it != polynomial.end())
        {
            // Monomial with the same exponents already exists, sum the coefficients
            it->second.setCoefficient(it->second.getCoefficient() + monomial.getCoefficient());
        }
        else
        {
            // Monomial with the given exponents does not exist, add it to the map
            polynomial[monomial.getExponents()] = monomial;
        }
    }

    // Overload * operator to compute the product of two polynomials
    Polynomial<Dimension> operator*(const Polynomial<Dimension> &other) const
    {
        Polynomial<Dimension> result;

        for (const auto &monomial1 : polynomial)
        {
            for (const auto &monomial2 : other.polynomial)
            {
                Monomial<Dimension> productMonomial = monomial1.second * monomial2.second;
                result.addMonomial(productMonomial);
            }
        }

        return result;
    }

    // Overload * operator to compute the product of a polynomial and a monomial
    Polynomial<Dimension> operator*(const Monomial<Dimension> &other) const
    {
        Polynomial<Dimension> result;

        for (const auto &monomial : polynomial)
        {
            Monomial<Dimension> productMonomial = monomial.second * other;
            result.addMonomial(productMonomial);
        }

        return result;
    }

    // Size
    size_t size() const
    {
        return polynomial.size();
    }

    // Getter
    const std::map<std::vector<unsigned int>, Monomial<Dimension>> &getPolynomial() const
    {
        return polynomial;
    }
};

/*
class QuadrinomialPower
{
private:
    // Vectors of monomials being the expansions of (ax+by+cz+d)^n
    std::vector<Monomial3D> quadrinomialExpanded;

public:
    QuadrinomialPower(const real &a, const real &b, const real &c, const real &d, const unsigned int &power)
    {
        for (unsigned int k1 = 0; k1 <= power; ++k1)
        {
            for (unsigned int k2 = 0; k2 <= power - k1; ++k2)
            {
                for (unsigned int k3 = 0; k3 <= power - k1 - k2; ++k3)
                {
                    unsigned int k4 = power - k1 - k2 - k3;
                    std::cout << k1 << " " << k2 << " " << k3 << " " << k4 << std::endl;
                    real coeff = multinomialCoeff(power, {k1, k2, k3, k4}) * std::pow(a, k1) * std::pow(b, k2) * std::pow(c, k3) * std::pow(d, k4);
                    quadrinomialExpanded.emplace_back(Monomial3D(k1, k2, k3, coeff));
                }
            }
        }
    }

    // Size
    size_t size() const
    {
        return quadrinomialExpanded.size();
    }

    // Getter
    const Monomial3D &getMonomial(const std::size_t &i) const
    {
        return quadrinomialExpanded[i];
    }
};
*/

class LinearTrinomialPower : public Polynomial<2>
{
public:
    LinearTrinomialPower(const real &a, const real &b, const real &c, const unsigned int &power)
    {
        for (unsigned int k1 = 0; k1 <= power; ++k1)
        {
            for (unsigned int k2 = 0; k2 <= power - k1; ++k2)
            {
                unsigned int k3 = power - k1 - k2;
                // std::cout << k1 << " " << k2 << " " << k3 << std::endl;
                real coeff = multinomialCoeff(power, {k1, k2, k3}) * std::pow(a, k1) * std::pow(b, k2) * std::pow(c, k3);
                this->addMonomial(Monomial2D(k1, k2, coeff));
            }
        }
    }
};

Polynomial<2> toPolynomial2D(const Monomial3D &m, const real &h, const real &h_F, const Point3D &X_P, const Point3D &X_F, const Point3D &e_x, const Point3D &e_y, const Point3D &e_z)
{
    /*
    std::cout << m << std::endl
              << X_P << std::endl
              << X_F << std::endl
              << h << std::endl
              << h_F << std::endl
              << e_x << std::endl
              << e_y << std::endl
              << e_z << std::endl;
    */
    real a, b, c;
    Point3D nX(Point3D(1, 0, 0).dot(e_x), Point3D(1, 0, 0).dot(e_y), Point3D(1, 0, 0).dot(e_z));
    a = h_F / h * nX[0];
    b = h_F / h * nX[1];
    c = (X_F - X_P)[0] / h;
    LinearTrinomialPower xi(a, b, c, m.getExponents()[0]);
    /*
    std::cout << "xi" << std::endl;
    for (const auto &monomialPair : xi.getPolynomial())
    {
        std::cout << monomialPair.second << std::endl;
    }
    */
    Point3D nY(Point3D(0, 1, 0).dot(e_x), Point3D(0, 1, 0).dot(e_y), Point3D(0, 1, 0).dot(e_z));
    a = h_F / h * nY[0];
    b = h_F / h * nY[1];
    c = (X_F - X_P)[1] / h;
    LinearTrinomialPower eta(a, b, c, m.getExponents()[1]);
    /*
    std::cout << "eta" << std::endl;
    for (const auto &monomialPair : eta.getPolynomial())
    {
        std::cout << monomialPair.second << std::endl;
    }
    */
    Point3D nZ(Point3D(0, 0, 1).dot(e_x), Point3D(0, 0, 1).dot(e_y), Point3D(0, 0, 1).dot(e_z));
    a = h_F / h * nZ[0];
    b = h_F / h * nZ[1];
    c = (X_F - X_P)[2] / h;
    LinearTrinomialPower zeta(a, b, c, m.getExponents()[2]);
    /*
    std::cout << "zeta" << std::endl;
    for (const auto &monomialPair : zeta.getPolynomial())
    {
        std::cout << monomialPair.second << std::endl;
    }
    std::cout << std::endl;
    */
    /*
    Polynomial<2> xietazeta = (xi * eta) * zeta;
    std::cout << "xietazeta" << std::endl;
    for (const auto &monomialPair : xietazeta.getPolynomial())
    {
        std::cout << monomialPair.second << std::endl;
    }
    std::cout << std::endl;
    */
    return (xi * eta) * zeta;
}

#endif // __MONOMIAL_HPP_