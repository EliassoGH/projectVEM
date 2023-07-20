#ifndef __POLYGOON_HPP_
#define __POLYGON_HPP_
#include "traits.hpp"
#include "edge.hpp"
#include <vector>
#include <memory>
#include <iostream>
#include <stdexcept>

namespace geometry
{
    template <typename EdgeType>
    class Polygon
    {
    private:
        std::vector<EdgeType> edges;

    public:
        // Constructor
        Polygon(const std::vector<EdgeType> &edges) : edges(edges) {}

        // Add an edge to the polygon
        void addEdge(const EdgeType &edge)
        {
            edges.push_back(edge);
        }

        // Get the number of edges in the polygon
        size_t size() const
        {
            return edges.size();
        }

        // Get a specific edge by index
        const EdgeType &getEdge(size_t index) const
        {
            if (index >= edges.size())
            {
                throw std::out_of_range("Invalid index for Polygon.");
            }
            return edges[index];
        }
    };
}

#endif // __EDGE_HPP_