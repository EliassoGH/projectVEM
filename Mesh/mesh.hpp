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

using namespace geometry;

template <typename PointType, typename EdgeType, typename PolygonType, typename PolyhedronType>
class Mesh
{
private:
    std::vector<PointType> points;
    std::vector<EdgeType> edges;
    std::vector<PolygonType> polygons;
    std::vector<PolyhedronType> polyhedra;

public:
    Mesh() = default;

    // Constructor
    Mesh(const std::vector<PointType> &_points, const std::vector<EdgeType> &_edges,
         const std::vector<PolygonType> &_polygons, const std::vector<PolyhedronType> &_polyhedra)
        : points(_points.begin(), _points.end()),
          edges(_edges.begin(), _edges.end()),
          polygons(_polygons.begin(), _polygons.end()),
          polyhedra(_polyhedra.begin(), _polyhedra.end()) {}

    // Constructor reading entities from Gmsh .geo file
    Mesh(const std::string &filename)
    {
        std::ifstream inputFile(filename);
        if (!inputFile.is_open())
        {
            throw std::runtime_error("Failed to open the file: " + filename);
        }

        // Temporary maps to handle the non ordered entries in the .geo file
        // I cannot resize the vectors and then use random access to insert entities since
        // not all entities can have a default constructor (e.g., edge, which requires two existing points)
        std::map<std::size_t, PointType> PointMap;
        std::map<std::size_t, EdgeType> LineMap;
        std::map<std::size_t, PolygonType> LineLoopMap;
        std::map<std::size_t, std::size_t> PlaneSurfaceMap;
        std::map<std::size_t, PolyhedronType> SurfaceLoopMap;

        /*
                std::size_t numPoints(0), numEdges(0);
                while (std::getline(inputFile, line))
                {
                    // Count the number of points in the .geo file
                    if (line.find("Point") == 0)
                    {
                        numPoints++;
                        continue;
                    }
                    // Count the number of edges in the .geo file
                    if (line.find("Line") == 0)
                    {
                        numEdges++;
                        continue;
                    }
                }
                points.resize(numPoints);
                edges.resize(numEdges);
                std::cout << "number of points reserved in points: " << points.capacity() << std::endl;

                // Clear the EOF flag and set the read position to the beginning of the file
                inputFile.clear();
                inputFile.seekg(0);
        */

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
                PointMap.insert(std::make_pair(id, point));
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

                EdgeType edge(PointMap.at(point1Id), PointMap.at(point2Id));
                edge.setId(id);
                LineMap.insert(std::make_pair(id, edge));
            }
            else if (type == "Line Loop")
            {
                // Read the Line Loop id
                std::size_t id;
                issline >> id;
                issline.ignore(3); // ignore the three characters of the string ")={"
                PolygonType polygon{};

                // Read the Lines
                long int signedEdgeId;
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
                    if (signedEdgeId < 0)
                    {
                        std::size_t edgeId = static_cast<std::size_t>(std::abs(signedEdgeId));
                        polygon.addEdge(LineMap.at(edgeId), 1);
                    }
                    else
                    {
                        std::size_t edgeId = static_cast<std::size_t>(signedEdgeId);
                        polygon.addEdge(LineMap.at(edgeId), 0);
                    }
                }
                LineLoopMap.insert(std::make_pair(id, polygon));
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
                PlaneSurfaceMap.insert(std::make_pair(id, LineLoopId));
                // Set the internal id of the Line Loop to the id of the Plane Surface
                LineLoopMap[LineLoopId].setId(id);
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
                    if (signedPlaneSurfaceId < 0)
                    {
                        std::size_t planeSurfaceId = static_cast<std::size_t>(std::abs(signedPlaneSurfaceId));
                        polyhedron.addPolygon(LineLoopMap.at(PlaneSurfaceMap.at(planeSurfaceId)),1);
                    }
                    else
                    {
                        std::size_t planeSurfaceId = static_cast<std::size_t>(std::abs(signedPlaneSurfaceId));
                        polyhedron.addPolygon(LineLoopMap.at(PlaneSurfaceMap.at(planeSurfaceId)),0);
                    }
                }
                polyhedron.setId(id);
                SurfaceLoopMap.insert(std::make_pair(id, polyhedron));
            }
        }

        // Assign values from the maps to the data structure vectors
        for (const auto &pair : PointMap)
        {
            points.push_back(pair.second);
        }
        for (const auto &pair : LineMap)
        {
            edges.push_back(pair.second);
        }
        for (const auto &pair : PlaneSurfaceMap)
        {
            polygons.push_back(LineLoopMap.at(pair.second));
        }
        for (const auto &pair : SurfaceLoopMap)
        {
            polyhedra.push_back(pair.second);
        }

        inputFile.close();
    }

    // Method to get the number of points
    std::size_t numPoints() const
    {
        return points.size();
    }
    // Method to get the number of edges
    std::size_t numEdges() const
    {
        return edges.size();
    }
    // Method to get the number of polygons
    std::size_t numPolygons() const
    {
        return polygons.size();
    }
    // Method to get the number of polyhedra
    std::size_t numPolyhedra() const
    {
        return polyhedra.size();
    }

    // Getter for a point
    const PointType &getPoint(std::size_t index) const
    {
        if (index >= points.size())
        {
            throw std::out_of_range("Invalid index for point.");
        }
        return points[index];
    }

    // Getter for an edge
    const EdgeType &getEdge(std::size_t index) const
    {
        if (index >= edges.size())
        {
            throw std::out_of_range("Invalid index for edge.");
        }
        return edges[index];
    }

    // Getter for a polygon
    const PolygonType &getPolygon(std::size_t index) const
    {
        if (index >= polygons.size())
        {
            throw std::out_of_range("Invalid index for polygon.");
        }
        return polygons[index];
    }

    // Getter for a polyhedron
    const PolyhedronType &getPolyhedron(std::size_t index) const
    {
        if (index >= polyhedra.size())
        {
            throw std::out_of_range("Invalid index for polyhedron.");
        }
        return polyhedra[index];
    }
};

#endif // __MESH_HPP_