#ifndef __POINT_HPP_
#define __POINT_HPP_
#include "traits.hpp"
#include <array>
#include <iostream>
#include <memory>
// #include <functional> // only needed for out-of-class piecewiseMultiply
#include <cmath>

namespace geometry
{
    template <typename... Args>
    class Point
    {
    private:
        static std::size_t lastId;
        std::size_t id;

    protected:
        std::array<real, sizeof...(Args)> coordinates;

    public:
        // Default constructor to initialize coordinates to 0
        Point() : coordinates{0}, id(lastId++) {}

        Point(Args... args) : coordinates{static_cast<real>(args)...}, id(lastId++) {}

        constexpr std::size_t getDimension() const
        {
            return coordinates.size();
        }

        // Set id
        void setId(std::size_t _id)
        {
            id=_id;
            lastId=_id+1;
        }

        // Get id
        const std::size_t getId() const
        {
            return id;
        }

        /*
                // Setter to modify all coordinates of a point
                template <typename... OtherArgs>
                void setCoordinates(OtherArgs... args)
                {
                    static_assert(sizeof...(Args) == sizeof...(OtherArgs), "Invalid number of arguments.");
                    std::size_t index = 0;
                    for (const auto coordinate : {args...})
                    {
                        coordinates[index] = coordinate;
                        index++;
                    }
                }
        */
        const std::array<real, sizeof...(Args)>
        getCoordinates() const
        {
            return coordinates;
        }

        /*
                real &operator[](std::size_t index)
                {
                    if (index >= sizeof...(Args))
                    {
                        throw std::out_of_range("Invalid dimension index.");
                    }
                    return coordinates[index];
                }
        */
        const real &operator[](std::size_t index) const
        {
            if (index >= sizeof...(Args))
            {
                throw std::out_of_range("Invalid dimension index.");
            }
            return coordinates[index];
        }

        // Overloaded * operator to compute the scalar multiplication of a point (point*scalar)
        auto operator*(const real &scalar) const
        {
            return multiplyByScalar(scalar, std::make_index_sequence<sizeof...(Args)>());
        }

        // Overloaded * operator to support scalar * point multiplication
        friend auto operator*(const real &scalar, const Point &point)
        {
            return point * scalar;
        }

        // scalar multiplication of a point
        template <size_t... Indices>
        auto multiplyByScalar(real scalar, std::index_sequence<Indices...>) const
        {
            return Point<Args...>((coordinates[Indices] * scalar)...);
        }

        // Overloaded / operator to compute the scalar division of a point (point/scalar)
        auto operator/(const real &scalar) const
        {
            return multiplyByScalar(1. / scalar, std::make_index_sequence<sizeof...(Args)>());
        }

        // Overloaded + operator to compute the sum of two points
        template <typename... OtherArgs>
        auto operator+(const Point<OtherArgs...> &other) const
        {
            static_assert(sizeof...(Args) == sizeof...(OtherArgs), "Invalid number of arguments.");
            return binaryOperation<OtherArgs...>(other, std::plus<>(), std::index_sequence_for<Args...>());
        }

        // Overloaded - operator to compute the difference of two points
        template <typename... OtherArgs>
        auto operator-(const Point<OtherArgs...> &other) const
        {
            static_assert(sizeof...(Args) == sizeof...(OtherArgs), "Invalid number of arguments.");
            return binaryOperation<OtherArgs...>(other, std::minus<>(), std::index_sequence_for<Args...>());
        }

        // Piecewise multiplication of two points
        template <typename... OtherArgs>
        auto piecewiseMultiply(const Point<OtherArgs...> &other) const
        {
            static_assert(sizeof...(Args) == sizeof...(OtherArgs), "Invalid number of arguments.");
            return binaryOperation<OtherArgs...>(other, std::multiplies<>(), std::index_sequence_for<Args...>());
        }

        // Binary operation in class
        template <typename... OtherArgs, size_t... Indices, typename Operation>
        auto binaryOperation(const Point<OtherArgs...> &other, Operation operation, std::index_sequence<Indices...>) const
        {
            return Point<OtherArgs...>(operation(coordinates[Indices], other[Indices])...);
        }

        // Dot product in the form point1.dot(point2)
        template <typename... OtherArgs>
        auto dot(const Point<OtherArgs...> &other) const
        {
            static_assert(sizeof...(Args) == sizeof...(OtherArgs), "Invalid number of arguments.");

            real sum = 0.0;
            sumOfPowers(this->piecewiseMultiply(other), 1, sum, std::index_sequence_for<Args...>());
            return sum;
        }

        // In-class cross product calculation (3D points only) in the form point1.cross(point2)
        template <typename... OtherArgs>
        auto cross(const Point<OtherArgs...> &other) const
        {
            static_assert(sizeof...(Args) == 3 && sizeof...(OtherArgs) == 3, "Cross product is only defined for 3D points.");
            return Point((*this)[1] * other[2] - (*this)[2] * other[1], (*this)[2] * other[0] - (*this)[0] * other[2], (*this)[0] * other[1] - (*this)[1] * other[0]);
        }

        // Euclidean distance
        template <typename... OtherArgs>
        auto distance(const Point<OtherArgs...> &other) const
        {
            static_assert(sizeof...(Args) == sizeof...(OtherArgs), "Invalid number of arguments.");

            auto diff = (*this) - other;
            return std::sqrt(diff.dot(diff));
        }

        // Define the comparison function based on edge Ids
        bool operator<(const Point<Args...> &other) const
        {
            if (*this == other)
            {
                return false;
            }
            else
            {
                return id < other.id;
            }
        }

        // Custom definition of operator== for Point
        template <typename... OtherArgs>
        bool operator==(const Point<OtherArgs...> &other) const
        {
            static_assert(sizeof...(Args) == sizeof...(OtherArgs), "Invalid number of arguments.");

            return coordinates == other.coordinates;
        }

        // Output stream operator to stream coordinates
        friend std::ostream &operator<<(std::ostream &os, const Point<Args...> &point)
        {
            os << "Point " << point.getId() << ": (";
            for (std::size_t i = 0; i < point.getDimension(); ++i)
            {
                if (i > 0)
                    os << " ";
                os << point.getCoordinates()[i];
            }
            os << ")";
            return os;
        }
    };

    /*
        // Function to compute piecewise multiplication between two points
        template <typename... Args, typename... OtherArgs>
        auto piecewiseMultiply(const Point<Args...> &p1, const Point<OtherArgs...> &p2)
        {
            static_assert(sizeof...(Args) == sizeof...(OtherArgs), "Invalid number of arguments.");

            std::function<real(real, real)> operation = [](real a, real b)
            { return a * b; };
            return binaryOperation(p1, p2, operation, std::index_sequence_for<Args...>());
        }
    */

    // Binary operation out of class
    template <typename... Args, typename... OtherArgs, size_t... Indices, typename Operation>
    auto binaryOperation(const Point<Args...> &p1, const Point<OtherArgs...> &p2, Operation operation, std::index_sequence<Indices...>)
    {
        return Point<Args...>(operation(p1.getCoordinates()[Indices], p2[Indices])...);
    }

    // Euclidean distance out of class, distance(point1,point2)
    template <typename... Args, typename... OtherArgs>
    auto distance(const Point<Args...> &p1, const Point<OtherArgs...> &p2)
    {
        return p1.distance(p2);
    }

    // Dot product out of class, dot(point1,point2)
    template <typename... Args, typename... OtherArgs>
    auto dot(const Point<Args...> &p1, const Point<OtherArgs...> &p2)
    {
        return p1.dot(p2);
    }

    // Cross product out of class, cross(point1,point2)
    template <typename... Args, typename... OtherArgs>
    auto cross(const Point<Args...> &p1, const Point<OtherArgs...> &p2)
    {
        return p1.cross(p2);
    }

    // Computes the sum of the powers
    template <typename... Args, size_t... Indices>
    void sumOfPowers(const Point<Args...> &point, const auto &p, real &sum, std::index_sequence<Indices...>)
    {
        ((sum += std::pow(point[Indices], p)), ...);
    }

    // Initialize lastId
    template <typename... Args>
    std::size_t Point<Args...>::lastId = 0;

    using Point3D = Point<real, real, real>;

}

#endif // __POINT_HPP_