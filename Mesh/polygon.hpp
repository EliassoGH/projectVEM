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
        static IndexType lastId;
        IndexType id;
        bool orientation;
        std::vector<std::reference_wrapper<const EdgeType>> edges;
        Point3D outwardNormalArea;
        real diameter = 0.0;
        // Point3D centroid;

    public:
        Polygon() : Polygon({}){};

        // Constructor taking individual edges
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

        // Set orientation
        void setOrientation(bool _orientation)
        {
            orientation = _orientation;
        }

        // Add an edge and its direction to the polygon
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

        // Set id
        void setId(IndexType _id)
        {
            id = _id;
            lastId = _id + 1;
        }

        // Get id
        const IndexType &getId() const
        {
            return id;
        }

        // Get the number of edges in the polygon
        std::size_t numEdges() const
        {
            return edges.size();
        }

        // Get an edge by index
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
                // return edges[edges.size() - index - 1].get();
            }
        }

        // Access edges through []
        const EdgeType &operator[](IndexType index) const
        {
            return getEdge(index);
        }

        // Get original Edge
        const EdgeType &getPositiveEdge(std::size_t index) const
        {
            if (this->getEdge(index).getId() > 0)
                return this->getEdge(index);
            if (orientation == false)
                return edges[index].get().getOtherHalfEdge();
            else
                return edges[edges.size() - index - 1].get();
        }

        // Check if the edges are stored consistently
        bool areEdgesConsistent() const
        {
            if (edges.size() <= 1)
            {
                return true; // A single edge or no edges are always consistent
            }

            size_t numEdges = edges.size();
            for (size_t i = 0; i < numEdges; ++i)
            {
                // std::cout<<edges[1].get()<<std::endl;
                // std::cout<<edges[2].get()<<std::endl;
                // std::cout<<edges[3].get()<<std::endl;
                // std::cout<<edges[0].get()<<std::endl;
                EdgeType currentEdge = edges[i].get();
                // std::cout<<currentEdge<<std::endl;
                // std::cout<<"qui1"<<std::endl;
                EdgeType nextEdge = edges[(i + 1) % numEdges].get();
                // std::cout<<nextEdge<<std::endl;
                // std::cout<<"qui2"<<std::endl;

                if (currentEdge[1] != nextEdge[0])
                {
                    return false; // Inconsistent edges found
                }
            }

            return true; // All edges are consistent
        }

        // Compute properties
        void computeProperties()
        {
            computeOutwardNormalArea();
            computeDiameter();
            // computeCentroid();
        }

        // Compute outward normal unit vector
        void computeOutwardNormalArea()
        {
            outwardNormalArea = Point3D();
            for (std::size_t e = 0; e < edges.size(); e++)
            {
                outwardNormalArea = outwardNormalArea + (0.5 * cross(this->getEdge(e)[0], this->getEdge(e)[1]));
            }
            /*
            real nx(0.0), ny(0.0), nz(0.0);
            for (std::size_t e = 0; e < edges.size(); e++)
            {
                nx += ((edges[e].get()[0][1] - edges[e].get()[1][1]) * (edges[e].get()[0][2] + edges[e].get()[1][2]));
                ny += ((edges[e].get()[0][2] - edges[e].get()[1][2]) * (edges[e].get()[0][0] + edges[e].get()[1][0]));
                nz += ((edges[e].get()[0][0] - edges[e].get()[1][0]) * (edges[e].get()[0][1] + edges[e].get()[1][1]));
            }
            outwardNormalArea = Point3D(nx, ny, nz)/2.0;
            */
        }

        // Get outward normal unit vector
        const Point3D getOutwardNormal() const
        {
            if (outwardNormalArea == Point3D())
            {
                throw std::logic_error("Outward normal area has not been computed yet.");
            }
            return outwardNormalArea.normalize();
        }

        // Get first local axis e_x
        const Point3D get_e_x() const
        {
            return this->getEdge(0).getDirection();
        }

        // Get second local axis e_y
        const Point3D get_e_y() const
        {
            return cross(this->getOutwardNormal(), this->get_e_x());
        }

        // Get area
        real getArea() const
        {
            if (outwardNormalArea == Point3D())
            {
                throw std::logic_error("Outward normal area has not been computed yet.");
            }
            return outwardNormalArea.norm();
        }

        // Compute diameter
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

        // Get diameter
        real getDiameter() const
        {
            if (diameter == 0.0)
            {
                throw std::logic_error("Diameter has not been computed yet.");
            }
            return diameter;
        }
        /*
                // Compute centroid
                void computeCentroid()
                {
                    real cx(0.0), cy(0.0), cz(0.0);
                    if (this->getOutwardNormal()[]==
                    for (std::size_t e = 0; e < edges.size(); e++)
                    {
                        cx += ((edges[e].get()[0][0] + edges[e].get()[1][0]) * ((edges[e].get()[0][0]*edges[e].get()[1][1]) - (edges[e].get()[1][0]*edges[e].get()[0][1])));
                        cy += ((edges[e].get()[0][1] + edges[e].get()[1][1]) * ((edges[e].get()[0][0]*edges[e].get()[1][1]) - (edges[e].get()[1][0]*edges[e].get()[0][1])));
                        cz += ((edges[e].get()[0][2] + edges[e].get()[1][2]) * ((edges[e].get()[0][0]*edges[e].get()[1][2]) - (edges[e].get()[1][0]*edges[e].get()[0][2])));
                        std::cout<<"cx"<<cx<<cy<<cz<<std::endl;
                    }
                    centroid = Point3D(cx,cy,cz)/(6*this->getArea());
                    std::cout<<centroid<<std::endl;
                }

                // Get centroid
                const Point3D &getCentroid() const
                {
                    if (centroid == Point3D())
                    {
                        throw std::logic_error("Centroid has not been computed yet.");
                    }
                    return centroid;
                }
        */
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
    IndexType Polygon<EdgeType>::lastId = 1;

    using Polygon3D = Polygon<Edge3D>;
}

#endif // __POLYGON_HPP_