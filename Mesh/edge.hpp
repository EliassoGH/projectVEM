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
    class EdgeBase
    {
    protected:
        static std::size_t lastId;
        std::size_t id;
        std::array<std::reference_wrapper<const PointType>, 2> points;
        real length;

    public:
        // Constructor with two points
        EdgeBase(const PointType &point1, const PointType &point2) : points({point1, point2}), id(lastId++)
        {
            // Cannot create an edge composed of the same point
            if (&point1 == &point2)
            {
                throw std::invalid_argument("Cannot create an edge with the same point.");
            }
        }

        // Set id
        void setId(std::size_t _id)
        {
            id = _id;
            lastId = _id + 1;
        }

        // Get id
        const std::size_t getId() const
        {
            return id;
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

        // Get length
        const real getLength() const
        {
            return length;
        }

        // Define the comparison function based on edge Ids
        bool operator<(const EdgeBase<PointType> &other) const
        {
            // if edges are made by the same points, one cannot be < than the other
            if (*this == other)
            {
                return false;
            }
            // if the edges do not have the same points, compare their ids
            else
            {
                return id < other.id;
            }
        }

        // Equality operator for edges
        bool operator==(const EdgeBase<PointType> &other) const
        {
            // Compare two edges by comparing the points they contain (regardless of order)
            return ((points[0].get() == (other.points[0].get()) && points[1].get() == (other.points[1].get())) ||
                    (points[0].get() == (other.points[1].get()) && points[1].get() == (other.points[0].get())));
        }
    };

    // Initialize lastId
    template <typename PointType>
    std::size_t EdgeBase<PointType>::lastId = 0;

    // Derived class for a half edge
    template <typename PointType>
    class Edge : public EdgeBase<PointType>
    {
    private:
        bool flipped;
        PointType direction;

    public:

        // Constructor with two points
        Edge(const PointType &point1, const PointType &point2, bool _flipped = false) : EdgeBase<PointType>(point1, point2), flipped(_flipped)
        {
            update();
        }

        // Update the properties of the edge.
        // When modifying one point, the edge autmatically sees the change as they are references
        // It is not the case for other data structures
        void update()
        {
            this->length = distance(this->points[0].get(), this->points[1].get());
            if (flipped == false)
                direction = (this->points[1].get() - this->points[0].get()) / this->length;
            else
                direction = (this->points[0].get() - this->points[1].get()) / this->length;
        }

        // Constant getter
        const PointType &operator[](size_t index) const
        {
            if (index >= 2)
            {
                throw std::out_of_range("Invalid index for Edge.");
            }
            if (flipped == false)
                return this->points[index].get();
            else
                return this->points[1 - index].get();
        }

        // Get direction
        const PointType getDirection() const
        {
            return direction;
        }

        // Stream output operator for the Edge class
        friend std::ostream &operator<<(std::ostream &os, const EdgeBase<PointType> &edge)
        {
            if (edge.flipped == false)
            {
                os << "Edge " << edge.getId() << ": (" << edge[0].getId() << ", " << edge[1].getId() << ")";
            }
            else
            {
                os << "Edge " << edge.getId() << ": (" << edge[1].getId() << ", " << edge[0].getId() << ")";
            }
            return os;
        }
    };

    using Edge3D = Edge<Point3D>;
}

#endif // __EDGE_HPP_