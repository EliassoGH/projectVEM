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
    /*
        std::vector<PointType> points;
        std::vector<std::pair<EdgeType, EdgeType>> edges;
        std::vector<std::pair<PolygonType, PolygonType>> polygons;
        std::vector<PolyhedronType> polyhedra;
    */
    // Temporary maps to handle the non ordered entries in the .geo file
    // I cannot resize the vectors and then use random access to insert entities since
    // not all entities can have a default constructor (e.g., edge, which requires two existing points)
    std::map<std::size_t, PointType> points;
    std::map<IndexType, EdgeType> edges;
    std::map<IndexType, PolygonType> polygons;
    std::map<std::size_t, std::size_t> PlaneSurfaceMap;
    std::map<std::size_t, PolyhedronType> polyhedra;

public:
    Mesh() = default;

    // Constructor reading entities from Gmsh .geo file
    Mesh(const std::string &filename)
    {
        std::ifstream inputFile(filename);
        if (!inputFile.is_open())
        {
            throw std::runtime_error("Failed to open the file: " + filename);
        }

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
                }
                PolygonType polygon2(polygon1);
                polygon2.setOrientation(true);
                polygon1.computeProperties();
                polygon2.computeProperties();
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
                polygonsTEMP.at(-LineLoopId).setId(-id);
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
                polyhedra.insert(std::make_pair(id, polyhedron));
            }
        }

        inputFile.close();
    }
    /*
        // Constructor
        Mesh(const std::vector<PointType> &_points, const std::vector<EdgeType> &_edges,
             const std::vector<PolygonType> &_polygons, const std::vector<PolyhedronType> &_polyhedra)
            : points(_points.begin(), _points.end()),
              edges(_edges.begin(), _edges.end()),
              polygons(_polygons.begin(), _polygons.end()),
              polyhedra(_polyhedra.begin(), _polyhedra.end()) {}
    */
    /*
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
         std::map<std::size_t, std::pair<EdgeType, EdgeType>> LineMap;
         std::map<std::size_t, std::pair<PolygonType, PolygonType>> LineLoopMap;
         std::map<std::size_t, std::size_t> PlaneSurfaceMap;
         std::map<std::size_t, PolyhedronType> SurfaceLoopMap;

                 // This is to count entities, it was used to resize the vectors
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

                 // Create and insert pair of edges
                 LineMap.insert(std::make_pair(id, createPairEdges(PointMap.at(point1Id), PointMap.at(point2Id), id)));
             }
             else if (type == "Line Loop")
             {
                 // Read the Line Loop id
                 std::size_t id;
                 issline >> id;
                 issline.ignore(3); // ignore the three characters of the string ")={"
                 PolygonType polygon1{};

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
                     if (signedEdgeId < 0)
                     {
                         std::size_t edgeId = static_cast<std::size_t>(std::abs(signedEdgeId));
                         polygon1.addEdge(LineMap.at(edgeId).second);
                     }
                     else
                     {
                         std::size_t edgeId = static_cast<std::size_t>(signedEdgeId);
                         polygon1.addEdge(LineMap.at(edgeId).first);
                     }
                 }
                 PolygonType polygon2(polygon1);
                 polygon2.setOrientation(true);
                 LineLoopMap.insert(std::make_pair(id, std::make_pair(polygon1, polygon2)));
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
                 LineLoopMap.at(LineLoopId).first.setId(id);
                 LineLoopMap.at(LineLoopId).second.setId(-id);
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
                     // If PlaneSurface Id is negative, in Gmsh cyclic edges of polygon correspond to inward normal to polyhedron
                     if (signedPlaneSurfaceId < 0)
                     {
                         std::size_t planeSurfaceId = static_cast<std::size_t>(std::abs(signedPlaneSurfaceId));
                         polyhedron.addPolygon(LineLoopMap.at(PlaneSurfaceMap.at(planeSurfaceId)).first);
                     }
                     else
                     {
                         std::size_t planeSurfaceId = static_cast<std::size_t>(signedPlaneSurfaceId);
                         polyhedron.addPolygon(LineLoopMap.at(PlaneSurfaceMap.at(planeSurfaceId)).second);
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
             // points.push_back(std::move(pair.second));
         }
         for (const auto &pair : LineMap)
         {
             edges.push_back(pair.second);
         }
         for (const auto &pair : PlaneSurfaceMap)
         {
             polygons.push_back(LineLoopMap.at(pair.second));
             // polygons.push_back(std::move(LineLoopMap.at(pair.second)));
         }
         for (const auto &pair : SurfaceLoopMap)
         {
             polyhedra.push_back(pair.second);
             // polyhedra.push_back(std::move(pair.second));
         }

         inputFile.close();

         // std::cout << polygons[1].first.getEdge(0)[0];
         std::cout << edges[1].first.getOtherHalfEdge();
         std::cout << polygons[1].second[1][1];
         std::cout << polygons[3].second[2][0];
         std::cout << polygons[2].first[1][0];
         std::cout << polyhedra[1][2][0];
     }
    */

    // Method to get the number of vertices
    std::size_t
    numVertices() const
    {
        return points.size();
    }
    // Method to get the number of edges
    std::size_t numEdges() const
    {
        return edges.size() / 2;
    }
    // Method to get the number of polygons
    std::size_t numPolygons() const
    {
        return polygons.size() / 2;
    }
    // Method to get the number of polyhedra
    std::size_t numPolyhedra() const
    {
        return polyhedra.size();
    }

    // Getter for a vertex
    const PointType &getVertex(std::size_t index) const
    {
        if (index >= points.size())
        {
            throw std::out_of_range("Invalid index for point.");
        }
        return points.at(index);
    }

    // Getter for an edge
    const EdgeType &getEdge(IndexType index) const
    {
        if (static_cast<std::size_t>(std::abs(index)) >= edges.size())
        {
            throw std::out_of_range("Invalid index for edge.");
        }
        return edges.at(index);
    }

    // Getter for a polygon
    const PolygonType &getPolygon(IndexType index) const
    {
        if (static_cast<std::size_t>(std::abs(index)) >= polygons.size())
        {
            throw std::out_of_range("Invalid index for polygon.");
        }
        return polygons.at(index);
    }

    // Getter for a polyhedron
    const PolyhedronType &getPolyhedron(std::size_t index) const
    {
        if (index >= polyhedra.size())
        {
            throw std::out_of_range("Invalid index for polyhedron.");
        }
        return polyhedra.at(index);
    }

    // Getter for the map of vertices
    const std::map<std::size_t, PointType> &getVertices() const
    {
        return points;
    }
/*
    // Getter for the map of edges
    const std::map<IndexType, EdgeType> &getEdges() const
    {
        return edges;
    }
*/
    // Getter for the map of edges
    const std::map<IndexType, EdgeType> getEdges() const
    {
        std::map<IndexType, EdgeType> positiveEdges;

        // Copy only the edges with positive indices to the new map
        std::copy_if(edges.begin(), edges.end(), std::inserter(positiveEdges, positiveEdges.end()),
                     [](const auto &edgePair)
                     { return edgePair.first > 0; });

        return positiveEdges;
    }
/*
    // Getter for the map of polygons
    const std::map<IndexType, PolygonType> &getPolygons() const
    {
        return polygons;
    }
*/
    // Getter for the map of polygons
    const std::map<IndexType, PolygonType> getPolygons() const
    {
        std::map<IndexType, PolygonType> positivePolygons;

        // Copy only the polygons with positive indices to the new map
        std::copy_if(polygons.begin(), polygons.end(), std::inserter(positivePolygons, positivePolygons.end()),
                     [](const auto &polygonPair)
                     { return polygonPair.first > 0; });

        return positivePolygons;
    }

    // Getter for the map of polyhedra
    const std::map<std::size_t, PolyhedronType> &getPolyhedra() const
    {
        return polyhedra;
    }
};

#endif // __MESH_HPP_