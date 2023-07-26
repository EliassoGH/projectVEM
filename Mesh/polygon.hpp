#ifndef __POLYGON_HPP_
#define __POLYGON_HPP_
#include "traits.hpp"
#include "edge.hpp"
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace geometry
{
    template <typename EdgeType>
    class Polygon
    {
    private:
        static std::size_t lastId;
        std::size_t id;
        std::vector<std::reference_wrapper<const EdgeType>> edges;
        std::vector<bool> edgeDirections; // Vector to store the reading direction of each edge

    public:
        Polygon() : id(lastId++) {};

        // Constructor taking individual edges
        Polygon(const std::initializer_list<EdgeType> &edgesWithoutDirection) : id(lastId++)
        {
            for (const auto &edge : edgesWithoutDirection)
            {
                for (const auto &existingEdge : edges)
                {
                    if (existingEdge.get() == edge)
                    {
                        throw std::invalid_argument("Cannot add coinciding edges to the Polygon.");
                    }
                }
                edges.push_back(std::cref(edge));
                edgeDirections.push_back(false); // Default direction is assumed to be false
            }
        }

        // Add an edge and its direction to the polygon
        void addEdge(const EdgeType &edge, const bool &direction = false)
        {
            for (const auto &existingEdge : edges)
            {
                if (existingEdge.get() == edge)
                {
                    throw std::invalid_argument("Cannot add coinciding edges to the Polygon.");
                }
            }
            edges.push_back(std::cref(edge));
            edgeDirections.push_back(direction);
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

        // Get the number of edges in the polygon
        std::size_t numEdges() const
        {
            return edges.size();
        }

        // Get an edge by index
        const EdgeType &getEdge(size_t index) const
        {
            if (index >= edges.size())
            {
                throw std::out_of_range("Invalid index for Polygon edge.");
            }
            return edges[index].get();
        }

        // Access edges through []
        const EdgeType &operator[](size_t index) const
        {
            return getEdge(index);
        }

        // Check if the edges are stored consistently
        bool areEdgesConsistent() const
        {
            if (edges.size() <= 1)
            {
                return true; // A single edge or no edges are always consistent
            }

            size_t numEdges = edges.size();
            size_t currentEdgeFlipped, nextEdgeFlipped;
            for (size_t i = 0; i < numEdges; ++i)
            {
                EdgeType currentEdge = edges[i].get();
                EdgeType nextEdge = edges[(i + 1) % numEdges].get();
                currentEdgeFlipped = 0;
                nextEdgeFlipped = 0;

                if (edgeDirections[i] == true)
                {
                    currentEdgeFlipped++;
                }
                if (edgeDirections[(i + 1) % numEdges] == true)
                {
                    nextEdgeFlipped++;
                }

                if (currentEdge[1 - currentEdgeFlipped] != nextEdge[0 + nextEdgeFlipped])
                {
                    return false; // Inconsistent edges found
                }
            }

            return true; // All edges are consistent
        }

        // Define the comparison function based on polygon Ids
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

        // Comparison operator for polygons (==)
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

        // Stream output operator for the Polygon class
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

    // Define the comparison function for reference_wrappers
    template <typename EdgeType>
    bool operator<(const std::reference_wrapper<const Polygon<EdgeType>> &lhs,
                   const std::reference_wrapper<const Polygon<EdgeType>> &rhs)
    {
        return lhs.get() < rhs.get();
    }

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
    std::size_t Polygon<EdgeType>::lastId = 0;

    using Polygon3D = Polygon<Edge3D>;
}

#endif // __POLYGON_HPP_