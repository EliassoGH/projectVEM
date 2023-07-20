#ifndef __EDGE_HPP_
#define __EDGE_HPP_
#include "traits.hpp"
#include "point.hpp"
#include <array>
#include <functional>
// #include <set> // used for existingEdges data structure
#include <iostream>
#include <stdexcept>

namespace geometry
{
    template <typename PointType>
    class Edge
    {
    private:
        std::array<std::reference_wrapper<const PointType>, 2> points;
        real length;
        PointType direction;

    public:
        // Constructor with two points
        Edge(const PointType &point1, const PointType &point2) : points({point1, point2})
        {
            // Cannot create an edge composed of the same point
            if (&point1 == &point2)
            {
                throw std::invalid_argument("Cannot create an edge with the same point.");
            }

            update();
        }

        void update()
        {
            length = distance(points[0].get(), points[1].get());
            direction = (points[1].get() - points[0].get()) / length;
        }
        
        /*
                // Getter
                PointType &operator[](size_t index)
                {
                    if (index >= 2)
                    {
                        throw std::out_of_range("Invalid index for Edge.");
                    }
                    return points[index];
                }
        */
        // Constant getter
        const PointType &operator[](size_t index) const
        {
            if (index >= 2)
            {
                throw std::out_of_range("Invalid index for Edge.");
            }
            return points[index].get();
        }

        // Get direction
        const PointType getDirection() const
        {
            return direction;
        }

        // Get length
        const real getLength() const
        {
            return length;
        }

        // Equality operator for edges
        template <typename OtherPointType>
        bool operator==(const Edge<OtherPointType> &other) const
        {
            // Compare two edges by comparing the points they contain (regardless of order)
            return ((points[0].get() == (other.points[0].get()) && points[1].get() == (other.points[1].get())) ||
                    (points[0].get() == (other.points[1].get()) && points[1].get() == (other.points[0].get())));
        }
    };

}

#endif // __EDGE_HPP_