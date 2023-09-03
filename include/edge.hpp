#ifndef __EDGE_HPP_
#define __EDGE_HPP_
#include "traits.hpp"
#include "point.hpp"
#include <array>
#include <functional>
#include <memory>
#include <iostream>
#include <stdexcept>

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
        std::shared_ptr<Edge<PointType>> otherHalfEdge; // Pointer to store the other half-edge

    public:
        /**
         * @brief Constructor with two points
         * 
         * @param point1 
         * @param point2 
         * @param _flipped if reading direction is flipped
         */
        Edge(const PointType &point1, const PointType &point2, bool _flipped = false) : id(lastId++), flipped(_flipped), points({point1, point2})
        {
            // Cannot create an edge composed of the same point
            if (&point1 == &point2)
            {
                throw std::invalid_argument("Cannot create an edge with the same point.");
            }

            update();
        }
        
        /**
         * @brief Setter method to set the other half-edge
         * 
         * @param otherEdge 
         */
        void setOtherHalfEdge(const Edge<PointType> &otherEdge)
        {
            otherHalfEdge = std::make_shared<Edge<PointType>>(otherEdge);
        }

        /**
         * @brief Getter method to get the other half-edge
         * 
         * @return Edge<PointType>& 
         */
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

        /**
         * @brief Update the properties of the edge.
                  When modifying one point, the edge autmatically
                  sees the change as they are references.
                  It is not the case for other data structures
         * 
         */
        void update()
        {
            length = distance(points[0].get(), points[1].get());
            if (flipped == false)
                direction = (points[1].get() - points[0].get()) / length;
            else
                direction = (points[0].get() - points[1].get()) / length;
        }

        /**
         * @brief Set id
         * 
         * @param _id 
         */
        void setId(IndexType _id)
        {
            id = _id;
            lastId++;
        }

        /**
         * @brief Get id
         * 
         * @return const IndexType& 
         */
        const IndexType &getId() const
        {
            return id;
        }

        /**
         * @brief Get length
         * 
         * @return const real& 
         */
        const real &getLength() const
        {
            return length;
        }

        /**
         * @brief Get direction
         * 
         * @return const PointType 
         */
        const PointType getDirection() const
        {
            return direction;
        }

        /**
         * @brief Constant getter
         * 
         * @param index 
         * @return const PointType& 
         */
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

        /**
         * @brief Define the comparison function based on edge Ids
         * 
         * @param other 
         * @return true 
         * @return false 
         */
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

        /**
         * @brief Equality operator for edges
         * 
         * @param other 
         * @return true 
         * @return false 
         */
        bool operator==(const Edge<PointType> &other) const
        {
            // Compare two edges by comparing the points they contain (regardless of order)
            return ((points[0].get() == (other.points[0].get()) && points[1].get() == (other.points[1].get())) ||
                    (points[0].get() == (other.points[1].get()) && points[1].get() == (other.points[0].get())));
        }

        /**
         * @brief Stream output operator for the Edge class
         * 
         * @param os 
         * @param edge 
         * @return std::ostream& 
         */
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

    /**
     * @brief Method to construct a pair of half-edges and set their id
     * 
     * @tparam PointType 
     * @param point1 
     * @param point2 
     * @param id 
     * @return std::pair<Edge<PointType>, Edge<PointType>> 
     */
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