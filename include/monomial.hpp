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
#include <array>

using namespace traits;
using namespace geometry;

template <unsigned int Dimension>
class Monomial
{
private:
    std::vector<unsigned int> exponents;
    real coefficient;

public:
    /**
     * @brief Default constructor
     */
    Monomial() : exponents(Dimension, 0), coefficient(0.0) {}

    /**
     * @brief Constructor
     * 
     * @param exponents 
     * @param coefficient 
     */
    Monomial(std::vector<unsigned int> exponents, real coefficient) : exponents(std::move(exponents)), coefficient(coefficient) {}

    /**
     * @brief Getter for the exponents
     * 
     * @return const std::vector<unsigned int>& 
     */
    const std::vector<unsigned int> &getExponents() const
    {
        return exponents;
    }

    /**
     * @brief Getter for the coefficient
     * 
     * @return real 
     */
    real getCoefficient() const
    {
        return coefficient;
    }

    /**
     * @brief Setter for the coefficient
     * 
     * @param coeff 
     */
    void setCoefficient(real coeff)
    {
        coefficient = coeff;
    }

    /**
     * @brief Compute the product of two monomials
     * 
     * @param other 
     * @return Monomial<Dimension> 
     */
    Monomial<Dimension> operator*(const Monomial<Dimension> &other) const
    {
        std::vector<unsigned int> newExponents(exponents.size(), 0);

        std::transform(exponents.begin(), exponents.end(), other.exponents.begin(), newExponents.begin(), std::plus<unsigned int>());

        return Monomial<Dimension>(std::move(newExponents), coefficient * other.coefficient);
    }

    /**
     * @brief Compute the derivative with respect to a variable
     * 
     * @param variableIndex 
     * @return Monomial<Dimension> 
     */
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

    /**
     * @brief Evaluate the monomial at a point
     * 
     * @tparam PointType 
     * @param point 
     * @return real 
     */
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

    /**
     * @brief Get the monomial order
     * 
     * @return unsigned int 
     */
    unsigned int getOrder() const
    {
        return std::accumulate(exponents.begin(), exponents.end(), 0);
    }

    /**
     * @brief Output stream operator
     * 
     * @param os 
     * @param monomial 
     * @return std::ostream& 
     */
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

    /**
     * @brief Default constructor
     * 
     */
    Monomial2D() : Monomial<2>() {}

    /**
     * @brief Constructor with exponents and coefficient
     * 
     * @param expX 
     * @param expY 
     * @param coeff 
     */
    Monomial2D(unsigned int expX, unsigned int expY, double coeff) : Monomial<2>({expX, expY}, coeff) {}

    /**
     * @brief Constructor for implicit conversion from Monomial<2> to Monomial2D
     * 
     * @param monomial 
     */
    Monomial2D(const Monomial<2> &monomial)
        : Monomial<2>(monomial.getExponents(), monomial.getCoefficient()) {}

    /**
     * @brief Method to compute the monomials up to a given degree
     * 
     * @param order_ 
     */
    static void computeMonomialsUpToOrder(unsigned int order_);

    /**
     * @brief Method to get the monomials up to a given degree
     * 
     * @param order_ 
     * @return const std::vector<Monomial2D> 
     */
    static const std::vector<Monomial2D> getMonomialsOrdered(unsigned int order_);

    /**
     * @brief Method to compute the laplacians up to a given degree
     * 
     * @param order_ 
     */
    static void computeLaplaciansToMonomialsOrdered(unsigned int order_);

    /**
     * @brief Method to get the laplacians up to a given degree
     * 
     * @param order_ 
     * @return const std::vector<std::pair<std::pair<real, std::size_t>,
     * std::pair<real, std::size_t>>> 
     */
    static const std::vector<std::pair<std::pair<real, std::size_t>,
                                       std::pair<real, std::size_t>>>
    getLaplaciansToMonomialsOrdered(unsigned int order_);

    /**
     * @brief Method to compute the product of two monomials and return a new Monomial2D instance
     * 
     * @param other 
     * @return Monomial2D 
     */
    Monomial2D operator*(const Monomial2D &other) const
    {
        return Monomial::operator*(other);
    }

    /**
     * @brief Method to compute the derivative with respect to a variable and return a new Monomial2D instance
     * 
     * @param variableIndex 
     * @return Monomial2D 
     */
    Monomial2D derivative(unsigned int variableIndex) const
    {
        return Monomial::derivative(variableIndex);
    }

    /**
     * @brief Compute the derivative with respect to x
     * 
     * @return Monomial2D 
     */
    Monomial2D dx() const
    {
        return this->derivative(0);
    }

    /**
     * @brief Compute the derivative with respect to y
     * 
     * @return Monomial2D 
     */
    Monomial2D dy() const
    {
        return this->derivative(1);
    }

    /**
     * @brief Overriding the output stream operator to print "x, y" for Monomial2D
     * 
     * @param os 
     * @param monomial 
     * @return std::ostream& 
     */
    friend std::ostream &operator<<(std::ostream &os, const Monomial2D &monomial);
};

class Monomial3D : public Monomial<3>
{
private:
    static unsigned int order;
    static std::vector<Monomial3D> monomials_ordered;
    static std::vector<std::array<std::pair<real, std::size_t>, 3>>
        gradients_to_monomials_ordered;
    static std::vector<std::array<std::pair<real, std::size_t>, 3>>
        laplacians_to_monomials_ordered;

public:
    /**
     * @brief Default constructor
     * 
     */
    Monomial3D() : Monomial<3>() {}

    /**
     * @brief Constructor with exponents and coefficient
     * 
     * @param expX 
     * @param expY 
     * @param expZ 
     * @param coeff 
     */
    Monomial3D(unsigned int expX, unsigned int expY, unsigned int expZ, double coeff) : Monomial<3>({expX, expY, expZ}, coeff) {}

    /**
     * @brief Constructor for implicit conversion from Monomial<3> to Monomial3D
     * 
     * @param monomial 
     */
    Monomial3D(const Monomial<3> &monomial)
        : Monomial<3>(monomial.getExponents(), monomial.getCoefficient()) {}

    /**
     * @brief Method to compute the ordered monomials up to a given degree
     * 
     * @param order_ 
     */
    static void computeMonomialsUpToOrder(unsigned int order_);

    /**
     * @brief Method to get the ordered monomials up to a given degree
     * 
     * @param order_ 
     * @return const std::vector<Monomial3D> 
     */
    static const std::vector<Monomial3D>
    getMonomialsOrdered(unsigned int order_);

    /**
     * @brief Method to compute the gradients of the ordered monomials up to a given degree
     * 
     * @param order_ 
     */
    static void computeGradientsToMonomialsOrdered(unsigned int order_);

    /**
     * @brief ethod to get the gradients of the ordered monomials up to a given degree
     * 
     * @param order_ 
     * @return const std::vector<std::array<std::pair<real, std::size_t>, 3>> 
     */
    static const std::vector<std::array<std::pair<real, std::size_t>, 3>>
    getGradientsToMonomialsOrdered(unsigned int order_);

    /**
     * @brief Method to compute the laplacians of the ordered monomials up to a given degree
     * 
     * @param order_ 
     */
    static void computeLaplaciansToMonomialsOrdered(unsigned int order_);

    /**
     * @brief Method to get the laplacians of the ordered monomials up to a given degree
     * 
     * @param order_ 
     * @return const std::vector<std::array<std::pair<real, std::size_t>, 3>> 
     */
    static const std::vector<std::array<std::pair<real, std::size_t>, 3>>
    getLaplaciansToMonomialsOrdered(unsigned int order_);

    /**
     * @brief Method to compute the product of two monomials and return a new Monomial2D instance
     * 
     * @param other 
     * @return Monomial3D 
     */
    Monomial3D operator*(const Monomial3D &other) const
    {
        return Monomial::operator*(other);
    }

    /**
     * @brief Method to compute the derivative with respect to a variable and return a new Monomial2D instance
     * 
     * @param variableIndex 
     * @return Monomial3D 
     */
    Monomial3D derivative(unsigned int variableIndex) const
    {
        return Monomial::derivative(variableIndex);
    }

    /**
     * @brief Compute the derivative with respect to x
     * 
     * @return Monomial3D 
     */
    Monomial3D dx() const
    {
        return this->derivative(0);
    }

    /**
     * @brief Compute the derivative with respect to y
     * 
     * @return Monomial3D 
     */
    Monomial3D dy() const
    {
        return this->derivative(1);
    }

    /**
     * @brief Compute the derivative with respect to z
     * 
     * @return Monomial3D 
     */
    Monomial3D dz() const
    {
        return this->derivative(2);
    }

    /**
     * @brief Overriding the output stream operator to print "x, y, z" for Monomial3D
     * 
     * @param os 
     * @param monomial 
     * @return std::ostream& 
     */
    friend std::ostream &operator<<(std::ostream &os, const Monomial3D &monomial);
};

/**
 * @brief Computes the factorial
 * 
 * @param n 
 * @return unsigned int 
 */
inline unsigned int factorial(unsigned int n)
{
    return std::tgamma(n + 1); // Factorial using tgamma
}

/**
 * @brief Computes the multinomial coefficient n!/(k1!k2!...)
 * 
 * @param n 
 * @param k_values 
 * @return unsigned int 
 */
inline unsigned int multinomialCoeff(unsigned int n, const std::vector<unsigned int> &k_values)
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
    std::map<std::vector<unsigned int>, Monomial<Dimension>> polynomial;
    unsigned int order = 0;

public:
    /**
     * @brief Default constructor
     * 
     */
    Polynomial() = default;

    /**
     * @brief Constructor
     * 
     * @param monomials 
     */
    Polynomial(const std::vector<Monomial<Dimension>> &monomials)
    {
        for (const Monomial<Dimension> &monomial : monomials)
        {
            addMonomial(monomial);
        }
    }

    /**
     * @brief Method to add a monomial to the polynomial
     * 
     * @param monomial 
     */
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
            if (monomial.getOrder() > order)
                order = monomial.getOrder();
        }
    }

    /**
     * @brief Get the polynomial order
     * 
     * @return unsigned int 
     */
    unsigned int getOrder() const
    {
        return order;
    }

    /**
     * @brief Overload * operator to compute the product of two polynomials
     * 
     * @param other 
     * @return Polynomial<Dimension> 
     */
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

    /**
     * @brief Overload * operator to compute the product of a polynomial and a monomial
     * 
     * @param other 
     * @return Polynomial<Dimension> 
     */
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

    /**
     * @brief get the number of monomials making up the polynomial
     * 
     * @return size_t 
     */
    size_t size() const
    {
        return polynomial.size();
    }

    /**
     * @brief Getter for the polynomial
     * 
     * @return const std::map<std::vector<unsigned int>, Monomial<Dimension>>& 
     */
    const std::map<std::vector<unsigned int>, Monomial<Dimension>> &getPolynomial() const
    {
        return polynomial;
    }
};

class LinearTrinomialPower : public Polynomial<2>
{
public:
    /**
     * @brief Construct a new Linear Trinomial Power object of the kind (ax+by+c)^power
     * 
     * @param a 
     * @param b 
     * @param c 
     * @param power 
     */
    LinearTrinomialPower(const real &a, const real &b, const real &c, const unsigned int &power);
};

/**
 * @brief Convert a Monomial3D defined in the polyhedron reference frame to a Polynomial<2> defined in on a face
 * 
 * @param m monomial to be converted
 * @param h diameter of the polyhedron
 * @param h_F diameter of the face
 * @param X_P origin of the polyhedron coordinate system
 * @param X_F origin of the face coordinate system
 * @param e_x first local face axis
 * @param e_y second local face axis
 * @param e_z third local face axis
 * @return Polynomial<2> 
 */
Polynomial<2> toPolynomial2D(const Monomial3D &m, const real &h, const real &h_F, const Point3D &X_P, const Point3D &X_F, const Point3D &e_x, const Point3D &e_y, const Point3D &e_z);

#endif // __MONOMIAL_HPP_