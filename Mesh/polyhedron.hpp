#ifndef __POLYHEDRON_HPP_
#define __POLYHEDRON_HPP_
#include "traits.hpp"
#include "polygon.hpp"
#include <set>
#include <functional>
#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace geometry
{
    template <typename PolygonType>
    class Polyhedron
    {
    private:
        static std::size_t lastId;
        std::size_t id;
        std::vector<std::reference_wrapper<const PolygonType>> polygons;
        real diameter = 0.0;

    public:
        // Constructor taking individual polygons
        Polyhedron(const std::initializer_list<PolygonType> &polygonsWithoutDirection) : id(lastId++)
        {
            for (const auto &polygon : polygonsWithoutDirection)
            {
                for (const auto &existingPolygon : polygons)
                {
                    if (existingPolygon.get() == polygon)
                    {
                        throw std::invalid_argument("Cannot add coinciding polygons to the Polyhedron.");
                    }
                }
                polygons.push_back(std::cref(polygon));
            }
            if (polygons.size() > 0)
                computeDiameter();
        }

        // Add a polygon and its direction to the polyhedron
        void addPolygon(const PolygonType &polygon)
        {
            for (const auto &existingPolygon : polygons)
            {
                if (existingPolygon.get() == polygon)
                {
                    throw std::invalid_argument("Cannot add coinciding polygons to the Polyhedron.");
                }
            }
            polygons.push_back(std::cref(polygon));
        }

        // Set id
        void setId(std::size_t _id)
        {
            id = _id;
            lastId = _id + 1;
        }

        // Get id
        const std::size_t &getId() const
        {
            return id;
        }

        // Method to get the number of polygons in the Polyhedron
        std::size_t numPolygons() const
        {
            return polygons.size();
        }

        // Method to get a polygon by index
        const PolygonType &getPolygon(std::size_t index) const
        {
            if (index >= polygons.size())
            {
                throw std::out_of_range("Invalid index for Polygon in Polyhedron.");
            }

            return polygons[index].get();
        }

        // Access operator [] to get a polygon by index
        const PolygonType &operator[](std::size_t index) const
        {
            return getPolygon(index);
        }

        // Compute diameter
        void computeDiameter()
        {
            diameter = 0.0;
            std::set<Point3D> vertices;
            for (std::size_t f = 0; f < polygons.size(); f++)
            {
                for (std::size_t e = 0; e < this->getPolygon(f).numEdges(); e++)
                {
                    vertices.insert(this->getPolygon(f)[e][0]);
                }
            }
            for (auto it1 = vertices.begin(); it1 != vertices.end(); ++it1)
            {
                for (auto it2 = std::next(it1); it2 != vertices.end(); ++it2)
                {
                    real dist = distance(*it1, *it2);
                    if (dist > diameter)
                        diameter = dist;
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

        // Stream output operator for the Polyhedron class
        friend std::ostream &operator<<(std::ostream &os, const Polyhedron<PolygonType> &polyhedron)
        {
            os << "Polyhedron " << polyhedron.getId() << ": ";
            for (std::size_t i = 0; i < polyhedron.numPolygons(); ++i)
            {
                if (i > 0)
                    os << " ";
                os << polyhedron.getPolygon(i).getId();
            }
            return os;
        }
    };

    // Initialize lastId
    template <typename PolygonType>
    std::size_t Polyhedron<PolygonType>::lastId = 0;

}

#endif // __POLYHEDRON_HPP_