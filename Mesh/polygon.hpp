#ifndef __POLYGON_HPP_
#define __POLYGON_HPP_
#include "traits.hpp"
#include "edge.hpp"
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <memory>

namespace geometry
{
    template <typename EdgeType>
    class Polygon
    {
    private:
        static IndexType lastId;
        IndexType id;
        bool orientation;
        std::vector<std::reference_wrapper<const EdgeType>> edges;
        std::shared_ptr<Polygon<EdgeType>> otherPolygon; // Pointer to store the other polygon
        Point3D outwardNormalArea;
        real diameter = 0.0;

    public:
        /**
         * @brief Construct a new Polygon object
         * 
         */
        Polygon() : Polygon({}){};

        /**
         * @brief Constructor taking individual edges
         * 
         * @param edges_ 
         * @param _orientation 
         */
        Polygon(const std::initializer_list<EdgeType> &edges_, bool _orientation = false) : id(lastId++), orientation(_orientation)
        {
            for (const auto &edge : edges_)
            {
                for (const auto &existingEdge : edges)
                {
                    if (existingEdge.get() == edge)
                    {
                        throw std::invalid_argument("Cannot add coinciding edges to the Polygon.");
                    }
                }
                edges.push_back(std::cref(edge));
            }
            if (edges.size() > 0)
                computeProperties();
        }

        /**
         * @brief Setter method to set the other polygon
         * 
         * @param otherPolygon_ 
         */
        void setOtherPolygon(const Polygon<EdgeType> &otherPolygon_)
        {
            otherPolygon = std::make_shared<Polygon<EdgeType>>(otherPolygon_);
        }

        /**
         * @brief Getter method to get the other polygon
         * 
         * @return Polygon<EdgeType>& 
         */
        Polygon<EdgeType> &getOtherPolygon() const
        {
            if (otherPolygon)
            {
                return *otherPolygon;
            }
            else
            {
                throw std::runtime_error("Twin polygon not set!");
            }
        }

        /**
         * @brief Set orientation
         * 
         * @param _orientation 
         */
        void setOrientation(bool _orientation)
        {
            orientation = _orientation;
        }

        /**
         * @brief Add an edge and its direction to the polygon
         * 
         * @param edge 
         */
        void addEdge(const EdgeType &edge)
        {
            for (const auto &existingEdge : edges)
            {
                if (existingEdge.get() == edge)
                {
                    throw std::invalid_argument("Cannot add coinciding edges to the Polygon.");
                }
            }
            edges.push_back(std::cref(edge));
        }

        /**
         * @brief Set id
         * 
         * @param _id 
         */
        void setId(IndexType _id)
        {
            id = _id;
            lastId = _id + 1;
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
         * @brief Get the number of edges in the polygon
         * 
         * @return std::size_t 
         */
        std::size_t numEdges() const
        {
            return edges.size();
        }

        /**
         * @brief Get an edge by index
         * 
         * @param index 
         * @return const EdgeType& 
         */
        const EdgeType &getEdge(std::size_t index) const
        {
            if (index >= edges.size())
            {
                throw std::out_of_range("Invalid index for Polygon edge.");
            }
            if (orientation == false)
                return edges[index].get();
            else
            {
                return edges[edges.size() - index - 1].get().getOtherHalfEdge();
            }
        }

        /**
         * @brief Access edges through []
         * 
         * @param index 
         * @return const EdgeType& 
         */
        const EdgeType &operator[](IndexType index) const
        {
            return getEdge(index);
        }

        /**
         * @brief Get original Edge
         * 
         * @param index 
         * @return const EdgeType& 
         */
        const EdgeType &getPositiveEdge(std::size_t index) const
        {
            if (this->getEdge(index).getId() > 0)
                return this->getEdge(index);
            if (orientation == false)
                return edges[index].get().getOtherHalfEdge();
            else
                return edges[edges.size() - index - 1].get();
        }

        /**
         * @brief Check if the edges are stored consistently
         * 
         * @return true 
         * @return false 
         */
        bool areEdgesConsistent() const
        {
            if (edges.size() <= 1)
            {
                return true; // A single edge or no edges are always consistent
            }

            size_t numEdges = edges.size();
            for (size_t i = 0; i < numEdges; ++i)
            {
                EdgeType currentEdge = edges[i].get();
                EdgeType nextEdge = edges[(i + 1) % numEdges].get();

                if (currentEdge[1] != nextEdge[0])
                {
                    return false; // Inconsistent edges found
                }
            }

            return true; // All edges are consistent
        }

        /**
         * @brief Compute properties
         * 
         */
        void computeProperties()
        {
            computeOutwardNormalArea();
            computeDiameter();
        }

        /**
         * @brief Compute outward normal unit vector
         * 
         */
        void computeOutwardNormalArea()
        {
            outwardNormalArea = Point3D();
            for (std::size_t e = 0; e < edges.size(); e++)
            {
                outwardNormalArea = outwardNormalArea + (0.5 * cross(this->getEdge(e)[0], this->getEdge(e)[1]));
            }
        }

        /**
         * @brief Get outward normal unit vector
         * 
         * @return const Point3D 
         */
        const Point3D getOutwardNormal() const
        {
            if (outwardNormalArea == Point3D())
            {
                throw std::logic_error("Outward normal area has not been computed yet.");
            }
            return outwardNormalArea.normalize();
        }

        /**
         * @brief Get first local axis e_x
         * 
         * @return const Point3D 
         */
        const Point3D get_e_x() const
        {
            return this->getEdge(0).getDirection();
        }

        /**
         * @brief Get second local axis e_y
         * 
         * @return const Point3D 
         */
        const Point3D get_e_y() const
        {
            return cross(this->getOutwardNormal(), this->get_e_x());
        }

        /**
         * @brief Get area
         * 
         * @return real 
         */
        real getArea() const
        {
            if (outwardNormalArea == Point3D())
            {
                throw std::logic_error("Outward normal area has not been computed yet.");
            }
            return outwardNormalArea.norm();
        }

        /**
         * @brief Compute diameter
         * 
         */
        void computeDiameter()
        {
            diameter = 0.0;
            for (std::size_t i = 0; i < (edges.size() - 1); i++)
            {
                for (std::size_t j = i + 1; j < edges.size(); j++)
                {
                    if (distance(this->getEdge(i)[0], this->getEdge(j)[0]) > diameter)
                        diameter = distance(this->getEdge(i)[0], this->getEdge(j)[0]);
                }
            }
        }

        /**
         * @brief Get diameter
         * 
         * @return real 
         */
        real getDiameter() const
        {
            if (diameter == 0.0)
            {
                throw std::logic_error("Diameter has not been computed yet.");
            }
            return diameter;
        }
        
        /**
         * @brief Define the comparison function based on polygon Ids
         * 
         * @param other 
         * @return true 
         * @return false 
         */
        bool operator<(const Polygon<EdgeType> &other) const
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

        /**
         * @brief Comparison operator for polygons (==)
         * 
         * @param other 
         * @return true 
         * @return false 
         */
        bool operator==(const Polygon<EdgeType> &other) const
        {

            // Check if the edges of the current polygon are a cyclic permutation of the edges in the other polygon
            if (!areEquivalentCyclicPermutations(dereferenceVector(edges), dereferenceVector(other.edges)))
            {
                return false;
            }

            // Otherwise the polygons are equivalent
            return true;
        }

        /**
         * @brief Stream output operator for the Polygon class
         * 
         * @param os 
         * @param polygon 
         * @return std::ostream& 
         */
        friend std::ostream &operator<<(std::ostream &os, const Polygon<EdgeType> &polygon)
        {
            os << "Polygon " << polygon.getId() << ": ";
            for (std::size_t i = 0; i < polygon.numEdges(); ++i)
            {
                if (i > 0)
                    os << " ";
                os << polygon.getEdge(i).getId();
            }
            return os;
        }
    };

    /**
     * @brief Define the comparison function for reference_wrappers
     * 
     * @tparam EdgeType 
     * @param lhs 
     * @param rhs 
     * @return true 
     * @return false 
     */
    template <typename EdgeType>
    bool operator<(const std::reference_wrapper<const Polygon<EdgeType>> &lhs,
                   const std::reference_wrapper<const Polygon<EdgeType>> &rhs)
    {
        return lhs.get() < rhs.get();
    }

    /**
     * @brief dereference a vector of reference_wrapper
     * 
     * @tparam T 
     * @param v 
     * @return std::vector<T> 
     */
    template <typename T>
    std::vector<T> dereferenceVector(const std::vector<std::reference_wrapper<const T>> &v)
    {
        std::vector<T> result;
        result.reserve(v.size());
        for (const auto &elem : v)
        {
            result.push_back(elem.get());
        }
        return result;
    }

    /**
     * @brief check if two vectors are equivalent cyclic permutations
     * 
     * @tparam T 
     * @param v1 
     * @param v2 
     * @return true 
     * @return false 
     */
    template <typename T>
    bool areEquivalentCyclicPermutations(const std::vector<T> &v1,
                                         const std::vector<T> &v2)
    {
        if (v1.size() != v2.size())
        {
            return false;
        }

        // Concatenate v1 with itself
        auto concatenated = v1;
        concatenated.insert(concatenated.end(), v1.begin(), v1.end());

        // Check if v2 is a contiguous subsequence of the concatenated vector
        auto it = std::search(concatenated.begin(), concatenated.end(), v2.begin(), v2.end());
        if (it != concatenated.end())
        {
            return true;
        }

        // Check if v2 is a contiguous subsequence of the reversed concatenated vector
        std::reverse(concatenated.begin(), concatenated.end());
        it = std::search(concatenated.begin(), concatenated.end(), v2.begin(), v2.end());
        return it != concatenated.end();
    }

    // Initialize lastId
    template <typename EdgeType>
    IndexType Polygon<EdgeType>::lastId = 1;

    using Polygon3D = Polygon<Edge3D>;
}

#endif // __POLYGON_HPP_