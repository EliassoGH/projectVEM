#ifndef __EDGE_HPP_
#define __EDGE_HPP_
#include "traits.hpp"
#include "point.hpp"
#include <array>
#include <memory>
//#include <set> // used for existingEdges data structure
#include <iostream>
#include <stdexcept>

namespace geometry
{
    template <typename PointType>
    class Edge
    {
    private:
        std::array<std::shared_ptr<PointType>, 2> points;
        real length;
        PointType direction;

        // Static set to keep track of existing pairs of points
        // static std::set<std::pair<std::shared_ptr<PointType>, std::shared_ptr<PointType>>> existingEdges;

    public:
        // Constructor with two points
        Edge(const PointType &point1, const PointType &point2)
        {
            // Cannot create an edge composed of the same point
            if (point1.getSharedPtr() == point2.getSharedPtr())
            {
                throw std::invalid_argument("Cannot create an edge with the same point.");
            }
            /*
                        // Sort the points to ensure consistent ordering in the set
                        std::shared_ptr<PointType> pmin = std::min(point1.getSharedPtr(), point2.getSharedPtr());
                        std::shared_ptr<PointType> pmax = std::max(point1.getSharedPtr(), point2.getSharedPtr());
                        std::cout<<"Address of p1: "<<point1.getSharedPtr()<<", address of p2: "<<point2.getSharedPtr()<<std::endl;
                        std::cout<<"Address of pmin: "<<pmin<<", address of pmax: "<<pmax<<std::endl;

                        // Cannot create an edge already existing
                        if (existingEdges.find(std::make_pair(pmin, pmax)) != existingEdges.end())
                        {
                            throw std::invalid_argument("Edge between the given points already exists.");
                        }
            */
            points[0] = point1.getSharedPtr();
            points[1] = point2.getSharedPtr();

            // Add the current edge to the static set when constructed
            // existingEdges.insert(std::make_pair(pmin, pmax));

            update();
        }

        void update()
        {
            length = distance(*points[0], *points[1]);
            direction = (*points[1] - *points[0]) / length;
        }
        /*
                // Output stream operator to stream addresses of existing edges
                static void printExistingEdges()
                {
                    std::cout<<"Existing edges addresses are: ";
                    for (const auto &edge : existingEdges)
                    {
                        std::cout << "( "<<edge.first <<", "<<edge.second<< "), ";
                    }
                    std::cout<<std::endl;
                }
        */
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
            return *(points[index]);
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
            return ((*points[0] == *(other.points[0]) && *points[1] == *(other.points[1])) ||
                    (*points[0] == *(other.points[1]) && *points[1] == *(other.points[0])));
        }
    };

    // Static set instantiation to keep track of existing pairs of points
    // template <typename PointType>
    // std::set<std::pair<std::shared_ptr<PointType>, std::shared_ptr<PointType>>> Edge<PointType>::existingEdges;

}

#endif // __EDGE_HPP_