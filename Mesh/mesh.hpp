#ifndef __MESH_HPP_
#define __MESH_HPP_
#include "traits.hpp"
#include "point.hpp"
#include "edge.hpp"
#include "polygon.hpp"
#include "polyhedron.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <utility>

using namespace geometry;

template <typename PointType, typename EdgeType, typename PolygonType, typename PolyhedronType>
class Mesh
{
private:
    std::map<std::size_t, PointType> points;
    std::map<IndexType, EdgeType> edges;
    std::map<IndexType, PolygonType> polygons;
    std::map<std::size_t, std::size_t> PlaneSurfaceMap;
    std::map<std::size_t, PolyhedronType> polyhedra;
    real h; // size of the mesh

public:
    /**
     * @brief Default constructor
     * 
     */
    Mesh() = default;

    /**
     * @brief Constructor reading entities from Gmsh .geo file
     * 
     * @param filename 
     */
    Mesh(const std::string &filename)
    {
        std::ifstream inputFile(filename);
        if (!inputFile.is_open())
        {
            throw std::runtime_error("Failed to open the file: " + filename);
        }

        h=0.0;

        std::map<IndexType, PolygonType> polygonsTEMP;

        std::string line;
        while (std::getline(inputFile, line))
        {
            std::istringstream issline(line);
            // Read the type
            std::string type;
            std::getline(issline, type, '(');

            // Check if it's a Point entity
            if (type == "Point")
            {
                // Read the Point id
                std::size_t id;
                issline >> id;
                issline.ignore(3); // ignore the three characters of the string ")={"

                // Read the coordinates
                real x, y, z;
                issline >> x;
                issline.ignore(); // ignore the ','
                issline >> y;
                issline.ignore(); // ignore the ','
                issline >> z;

                PointType point(x, y, z);
                point.setId(id);
                points.insert(std::make_pair(id, point));
            }
            else if (type == "Line")
            {
                // Read the Line id
                std::size_t id;
                issline >> id;
                issline.ignore(3); // ignore the three characters of the string ")={"

                // Read the Points
                std::size_t point1Id, point2Id;
                issline >> point1Id;
                issline.ignore(); // ignore the ','
                issline >> point2Id;

                // Create and insert twin edges
                Edge<PointType> edge1(points.at(point1Id), points.at(point2Id), false);
                Edge<PointType> edge2(points.at(point1Id), points.at(point2Id), true);
                edge1.setId(id);
                edge2.setId(-id);
                edge1.setOtherHalfEdge(edge2);
                edge2.setOtherHalfEdge(edge1);
                edges.insert(std::make_pair(id, edge1));
                edges.insert(std::make_pair(-id, edge2));
            }
            else if (type == "Line Loop")
            {
                // Read the Line Loop id
                std::size_t id;
                issline >> id;
                issline.ignore(3); // ignore the three characters of the string ")={"
                PolygonType polygon1{};
                PolygonType polygon2{};

                // Read the Lines
                IndexType signedEdgeId;
                std::string edgeIdStr;
                while (std::getline(issline, edgeIdStr, ','))
                {
                    // Remove part coming from '}' on
                    std::size_t bracketFound = edgeIdStr.find('}');
                    if (bracketFound != std::string::npos)
                    {
                        edgeIdStr = edgeIdStr.substr(0, bracketFound);
                    }
                    signedEdgeId = std::stoll(edgeIdStr);
                    polygon1.addEdge(edges.at(signedEdgeId));
                    polygon2.addEdge(edges.at(signedEdgeId));
                }
                //PolygonType polygon2(polygon1);
                polygon2.setOrientation(true);
                polygon1.computeProperties();
                polygon2.computeProperties();
                polygon1.setOtherPolygon(polygon2);
                polygon2.setOtherPolygon(polygon1);
                polygonsTEMP.insert(std::make_pair(id, polygon1));
                polygonsTEMP.insert(std::make_pair(-id, polygon2));
            }
            else if (type == "Plane Surface")
            {
                // Read the Plane Surface id, which will be the Polygon id
                std::size_t id;
                issline >> id;
                issline.ignore(3); // ignore the three characters of the string ")={"

                // Read the Line Loop
                std::size_t LineLoopId;
                issline >> LineLoopId;

                polygonsTEMP.at(LineLoopId).setId(id);
                polygonsTEMP.at(LineLoopId).getOtherPolygon().setId(-id);
                polygonsTEMP.at(-LineLoopId).setId(-id);
                polygonsTEMP.at(-LineLoopId).getOtherPolygon().setId(id);
                polygons.insert(std::make_pair(id, polygonsTEMP.at(LineLoopId)));
                polygons.insert(std::make_pair(-id, polygonsTEMP.at(-LineLoopId)));
            }
            else if (type == "Surface Loop")
            {
                // Read the Surface Loop id
                std::size_t id;
                issline >> id;
                issline.ignore(3); // ignore the three characters of the string ")={"
                PolyhedronType polyhedron{};

                // Read the Plane Surfaces
                long int signedPlaneSurfaceId;
                std::string planeSurfaceIdStr;
                while (std::getline(issline, planeSurfaceIdStr, ','))
                {
                    // Remove part coming from '}' on
                    std::size_t bracketFound = planeSurfaceIdStr.find('}');
                    if (bracketFound != std::string::npos)
                    {
                        planeSurfaceIdStr = planeSurfaceIdStr.substr(0, bracketFound);
                    }
                    signedPlaneSurfaceId = std::stoll(planeSurfaceIdStr);
                    polyhedron.addPolygon(polygons.at(-signedPlaneSurfaceId));
                }
                polyhedron.setId(id);
                polyhedron.computeDiameter();
                h+=polyhedron.getDiameter();
                polyhedra.insert(std::make_pair(id, polyhedron));
            }
        }
        h/=polyhedra.size();

        inputFile.close();
    }

    /**
     * @brief Method to get the number of vertices
     * 
     * @return std::size_t 
     */
    std::size_t
    numVertices() const
    {
        return points.size();
    }

    /**
     * @brief Method to get the number of edges
     * 
     * @return std::size_t 
     */
    std::size_t numEdges() const
    {
        return edges.size() / 2;
    }

    /**
     * @brief Method to get the number of polygons
     * 
     * @return std::size_t 
     */
    std::size_t numPolygons() const
    {
        return polygons.size() / 2;
    }

    /**
     * @brief Method to get the number of polyhedra
     * 
     * @return std::size_t 
     */
    std::size_t numPolyhedra() const
    {
        return polyhedra.size();
    }

    /**
     * @brief Getter for a vertex
     * 
     * @param index 
     * @return const PointType& 
     */
    const PointType &getVertex(std::size_t index) const
    {
        if (index >= points.size())
        {
            throw std::out_of_range("Invalid index for point.");
        }
        return points.at(index);
    }

    /**
     * @brief Getter for an edge
     * 
     * @param index 
     * @return const EdgeType& 
     */
    const EdgeType &getEdge(IndexType index) const
    {
        if (static_cast<std::size_t>(std::abs(index)) >= edges.size())
        {
            throw std::out_of_range("Invalid index for edge.");
        }
        return edges.at(index);
    }

    /**
     * @brief Getter for a polygon
     * 
     * @param index 
     * @return const PolygonType& 
     */
    const PolygonType &getPolygon(IndexType index) const
    {
        if (static_cast<std::size_t>(std::abs(index)) >= polygons.size())
        {
            throw std::out_of_range("Invalid index for polygon.");
        }
        return polygons.at(index);
    }

    /**
     * @brief Getter for a polyhedron
     * 
     * @param index 
     * @return const PolyhedronType& 
     */
    const PolyhedronType &getPolyhedron(std::size_t index) const
    {
        if (index >= polyhedra.size())
        {
            throw std::out_of_range("Invalid index for polyhedron.");
        }
        return polyhedra.at(index);
    }

    /**
     * @brief Getter for the map of vertices
     * 
     * @return const std::map<std::size_t, PointType>& 
     */
    const std::map<std::size_t, PointType> &getVertices() const
    {
        return points;
    }
    
    /**
     * @brief Getter for the map of edges
     * 
     * @return const std::map<IndexType, EdgeType> 
     */
    const std::map<IndexType, EdgeType> getEdges() const
    {
        std::map<IndexType, EdgeType> positiveEdges;

        // Copy only the edges with positive indices to the new map
        std::copy_if(edges.begin(), edges.end(), std::inserter(positiveEdges, positiveEdges.end()),
                     [](const auto &edgePair)
                     { return edgePair.first > 0; });

        return positiveEdges;
    }

    /**
     * @brief Getter for the map of polygons
     * 
     * @return const std::map<IndexType, PolygonType> 
     */
    const std::map<IndexType, PolygonType> getPolygons() const
    {
        std::map<IndexType, PolygonType> positivePolygons;

        // Copy only the polygons with positive indices to the new map
        std::copy_if(polygons.begin(), polygons.end(), std::inserter(positivePolygons, positivePolygons.end()),
                     [](const auto &polygonPair)
                     { return polygonPair.first > 0; });

        return positivePolygons;
    }

    /**
     * @brief Getter for the map of polyhedra
     * 
     * @return const std::map<std::size_t, PolyhedronType>& 
     */
    const std::map<std::size_t, PolyhedronType> &getPolyhedra() const
    {
        return polyhedra;
    }

    /**
     * @brief Getter for the average diameter
     * 
     * @return const real& 
     */
    const real &getSize() const
    {
        return h;
    }

    /**
     * @brief Print mesh information
     * 
     */
    void print() const
    {
        std::cout<<"MESH"<<std::endl;
        std::cout<<"Average element size            : "<<this->getSize()<<std::endl;
        std::cout<<"Number of vertices              : "<<this->numVertices()<<std::endl;
        std::cout<<"Number of edges                 : "<<this->numEdges()<<std::endl;
        std::cout<<"Number of polygons              : "<<this->numPolygons()<<std::endl;
        std::cout<<"Number of polyhedra             : "<<this->numPolyhedra()<<std::endl;
    }
};

#endif // __MESH_HPP_