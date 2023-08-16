#ifndef __EDGE_HPP_
#define __EDGE_HPP_
#include "traits.hpp"
#include "point.hpp"
#include <array>
#include <functional>
#include <memory>
// #include <set> // used for existingEdges data structure
#include <iostream>
#include <stdexcept>
#include <optional>

namespace geometry
{
    template <typename PointType>
    class Edge
    {
    private:
        static IndexType lastId; // current available id
        IndexType id;
        bool flipped;
        std::array<std::reference_wrapper<const PointType>, 2> points;
        real length;
        PointType direction;
        // Edge<PointType> *otherHalfEdge; // Pointer to store the other half-edge
        std::shared_ptr<Edge<PointType>> otherHalfEdge; // Pointer to store the other half-edge
        // std::weak_ptr<Edge<PointType>> otherHalfEdge; // Use std::weak_ptr for non-owning reference
        // std::optional<std::reference_wrapper<const Edge<PointType>>> otherHalfEdge;

    public:
        // Constructor with two points
        Edge(const PointType &point1, const PointType &point2, bool _flipped = false) : id(lastId++), flipped(_flipped), points({point1, point2})
        {
            // Cannot create an edge composed of the same point
            if (&point1 == &point2)
            {
                throw std::invalid_argument("Cannot create an edge with the same point.");
            }

            update();
        }
        /*
                // Method to set the other half-edge
                void setOtherHalfEdge(Edge<PointType> &otherEdge)
                {
                    otherHalfEdge = std::make_shared<Edge<PointType>>(otherEdge);
                }

                // Method to get the other half-edge
                std::shared_ptr<Edge<PointType>> getOtherHalfEdge() const
                // Edge<PointType> getOtherHalfEdge() const
                {
                    if (!otherHalfEdge)
                    {
                        throw std::runtime_error("Twin edge not set!");
                    }
                    return otherHalfEdge;
                }
        */
        /*
        void setOtherHalfEdge(const Edge<PointType> &otherEdge)
        {
            otherHalfEdge = std::ref(otherEdge);
            std::cout << "other half edge has just been set to" << std::endl;
            std::cout << otherHalfEdge.value() << std::endl;
        }

        const Edge<PointType> &getOtherHalfEdge() const
        {
            if (otherHalfEdge.has_value())
            {
                std::cout << std::endl<< "other half edge is set" << std::endl;
                // Return the reference to the other edge inside the optional
                return otherHalfEdge.value().get();
            }
            else
            {
                std::cout << "other half edge not set" << std::endl;
                // Return a reference to the edge itself
                return *this;
            }
        }

        void setOtherHalfEdge(const Edge<PointType> &otherEdge)
        {
            *otherHalfEdge = otherEdge;
        }

        const Edge<PointType> &getOtherHalfEdge() const
        {
            return *otherHalfEdge;
        }
*/
        // Setter method to set the other half-edge
        void setOtherHalfEdge(const Edge<PointType> &otherEdge)
        {
            otherHalfEdge = std::make_shared<Edge<PointType>>(otherEdge);
        }

        // Getter method to get the other half-edge
        Edge<PointType> &getOtherHalfEdge() const
        {
            if (otherHalfEdge)
            {
                return *otherHalfEdge;
            }
            else
            {
                throw std::runtime_error("Twin edge not set!");
            }
        }
        // Update the properties of the edge.
        // When modifying one point, the edge autmatically sees the change as they are references
        // It is not the case for other data structures
        void update()
        {
            length = distance(points[0].get(), points[1].get());
            if (flipped == false)
                direction = (points[1].get() - points[0].get()) / length;
            else
                direction = (points[0].get() - points[1].get()) / length;
        }

        // Set id
        void setId(IndexType _id)
        {
            id = _id;
            lastId++;
        }

        // Get id
        const IndexType &getId() const
        {
            return id;
        }

        // Get length
        const real &getLength() const
        {
            return length;
        }

        // Get direction
        const PointType getDirection() const
        {
            return direction;
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
        const PointType &operator[](IndexType index) const
        {
            if (index >= 2)
            {
                throw std::out_of_range("Invalid index for Edge.");
            }
            if (flipped == false)
                return points[index].get();
            else
                return points[1 - index].get();
        }

        // Define the comparison function based on edge Ids
        bool operator<(const Edge<PointType> &other) const
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
        bool operator==(const Edge<PointType> &other) const
        {
            // Compare two edges by comparing the points they contain (regardless of order)
            return ((points[0].get() == (other.points[0].get()) && points[1].get() == (other.points[1].get())) ||
                    (points[0].get() == (other.points[1].get()) && points[1].get() == (other.points[0].get())));
        }

        // Stream output operator for the Edge class
        friend std::ostream &operator<<(std::ostream &os, const Edge<PointType> &edge)
        {
            os << "Edge " << edge.getId() << ": (" << edge[0].getId() << ", " << edge[1].getId() << ")";
            return os;
        }
    };

    // Initialize lastId
    template <typename PointType>
    IndexType Edge<PointType>::lastId = 1;

    using Edge3D = Edge<Point3D>;

    // Method to construct a pair of half-edges and set their id
    template <typename PointType>
    std::pair<Edge<PointType>, Edge<PointType>> createPairEdges(const PointType &point1, const PointType &point2, IndexType id)
    {
        Edge<PointType> edge1(point1, point2, false);
        Edge<PointType> edge2(point1, point2, true);
        edge1.setId(id);
        edge2.setId(-id);
        edge1.setOtherHalfEdge(edge2);
        edge2.setOtherHalfEdge(edge1);
        return std::make_pair(edge1, edge2);
    }
}

#endif // __EDGE_HPP_