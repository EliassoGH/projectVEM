#ifndef __VIRTUALDOFS_HPP_
#define __VIRTUALDOFS_HPP_

#include "monomial.hpp"
#include "integration.hpp"
#include "mesh.hpp"
#include <memory>
#include <vector>
#include <utility>
#include <typeinfo>

// Forward declaration
class VertexDof;
class EdgeDof;
class FaceDof;
class PolyhedronDof;

// Base class for degrees of freedom
class VirtualDof
{
public:
    virtual ~VirtualDof() = default;
    virtual std::size_t getId() const = 0;
    virtual std::ostream &operator<<(std::ostream &os) const = 0;

    friend std::ostream &operator<<(std::ostream &os, const VertexDof &vDof);
    friend std::ostream &operator<<(std::ostream &os, const EdgeDof &eDof);
    friend std::ostream &operator<<(std::ostream &os, const FaceDof &fDof);
    friend std::ostream &operator<<(std::ostream &os, const PolyhedronDof &pDof);
};

// Derived class for VertexDof
class VertexDof : public VirtualDof
{
private:
    std::size_t vertexId;
    Point3D vertex;

public:
    VertexDof(std::size_t id, const Point3D &vertex_) : vertexId(id), vertex(vertex_) {}

    std::size_t getId() const override
    {
        return vertexId;
    }

    const Point3D &getVertex() const
    {
        return vertex;
    }

    std::ostream &operator<<(std::ostream &os) const override
    {
        os << "Vertex Dof - Vertex ID: " << getId() << ", Coordinates: " << getVertex() << "\n";
        return os;
    }
};

// Define the stream output operator outside the class definition
std::ostream &operator<<(std::ostream &os, const VertexDof &vDof)
{
    os << "Vertex Dof - Vertex ID: " << vDof.getId() << ", Coordinates: " << vDof.getVertex() << "\n";
    return os;
}

// Derived class for EdgeDof
class EdgeDof : public VirtualDof
{
private:
    std::size_t edgeId;
    Point3D gaussLobattoPoint;
    real weight;

public:
    EdgeDof(std::size_t id, const Point3D &gaussLobattoPoint_, const real &weight_)
        : edgeId(id), gaussLobattoPoint(gaussLobattoPoint_), weight(weight_) {}

    std::size_t getId() const override
    {
        return edgeId;
    }

    const Point3D &getGaussLobattoPoint() const
    {
        return gaussLobattoPoint;
    }

    const real &getWeight() const
    {
        return weight;
    }

    std::ostream &operator<<(std::ostream &os) const override
    {
        os << "Edge Dof - Edge ID: " << getId() << ", Gauss-Lobatto Point: " << getGaussLobattoPoint() << ", weight: " << getWeight() << "\n";
        return os;
    }
};

// Define the stream output operator outside the class definition
std::ostream &operator<<(std::ostream &os, const EdgeDof &eDof)
{
    os << "Edge Dof - Edge ID: " << eDof.getId() << ", Coordinates: " << eDof.getGaussLobattoPoint() << ", weight: " << eDof.getWeight() << "\n";
    return os;
}

// Derived class for FaceDof
class FaceDof : public VirtualDof
{
private:
    std::size_t faceId;
    Monomial2D monomial;

public:
    FaceDof(std::size_t id, const Monomial2D &monomial_) : faceId(id), monomial(monomial_) {}

    std::size_t getId() const override
    {
        return faceId;
    }

    const Monomial2D &getMonomial() const
    {
        return monomial;
    }

    std::ostream &operator<<(std::ostream &os) const override
    {
        os << "Face Dof - Face ID: " << getId() << ", Monomial: " << getMonomial() << "\n";
        return os;
    }
};

// Define the stream output operator outside the class definition
std::ostream &operator<<(std::ostream &os, const FaceDof &fDof)
{
    os << "Face Dof - Face ID: " << fDof.getId() << ", Monomial: " << fDof.getMonomial() << "\n";
    return os;
}

// Derived class for PolyhedronDof
class PolyhedronDof : public VirtualDof
{
private:
    std::size_t polyhedronId;
    Monomial3D monomial;

public:
    PolyhedronDof(std::size_t id, const Monomial3D &monomial_) : polyhedronId(id), monomial(monomial_) {}

    std::size_t getId() const override
    {
        return polyhedronId;
    }

    const Monomial3D &getMonomial() const
    {
        return monomial;
    }

    std::ostream &operator<<(std::ostream &os) const override
    {
        os << "Polyhedron Dof - Polyhedron ID: " << getId() << ", Monomial: " << getMonomial() << "\n";
        return os;
    }
};

// Define the stream output operator outside the class definition
std::ostream &operator<<(std::ostream &os, const PolyhedronDof &pDof)
{
    os << "Polyhedron Dof - Polyhedron ID: " << pDof.getId() << ", Monomial: " << pDof.getMonomial() << "\n";
    return os;
}

class VirtualDofsCollection
{
private:
    static unsigned int order;
    std::vector<std::shared_ptr<VirtualDof>> dofs;
    std::size_t numVdofs = 0;
    std::size_t numEdofs = 0;
    std::size_t numFdofs = 0;
    std::size_t numPdofs = 0;

public:
    virtual ~VirtualDofsCollection() = default;

    // Constructor to create VirtualDofs from a Mesh and order
    VirtualDofsCollection(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh, unsigned int order_)
    {
        order = order_;

        // Create vertexDof objects
        numVdofs = mesh.numVertices();
        for (const auto &vertexPair : mesh.getVertices())
        {
            std::size_t vertexId = vertexPair.first;
            const Point3D &vertex = vertexPair.second;

            // Create a new VertexDof object and add it to the dofs vector
            std::shared_ptr<VertexDof> vDof = std::make_shared<VertexDof>(vertexId, vertex);
            dofs.push_back(vDof);
        }

        if (order > 1)
        {
            // Create EdgeDof objects
            numEdofs = mesh.numEdges() * order;
            for (const auto &edgePair : mesh.getEdges())
            {
                std::size_t edgeId = edgePair.first;
                const Edge3D &edge = edgePair.second;

                // Compute Gauss-Lobatto points for the edge based on the order
                auto pointsWeights = GaussLobatto::computeGaussLobattoPointsWOnEdge(edge, order + 1);

                // Create EdgeDof instances for each Gauss-Lobatto internal point on the edge
                for (std::size_t i = 1; i < (pointsWeights.first.size() - 1); i++)
                {
                    std::shared_ptr<EdgeDof> eDof = std::make_shared<EdgeDof>(edgeId, pointsWeights.first[i], pointsWeights.second[i]);
                    dofs.push_back(eDof);
                }
            }

            // Create FaceDof objects
            numFdofs = mesh.numPolygons() * ((Monomial2D::getMonomialsOrdered(order - 2)).size());
            for (const auto &facePair : mesh.getPolygons())
            {
                std::size_t faceId = facePair.first;
                const Polygon3D &face = facePair.second;

                // Get the monomials in 2D up to order k-2, where k order of the VEM
                auto monomials2D = Monomial2D::getMonomialsOrdered(order - 2);

                // Create FaceDof instance for each dof
                for (const auto &monomial : monomials2D)
                {
                    std::shared_ptr<FaceDof> fDof = std::make_shared<FaceDof>(faceId, monomial);
                    dofs.push_back(fDof);
                }
            }

            // Create PolyhedronDof objects
            numPdofs = mesh.numPolyhedra() * ((Monomial3D::getMonomialsOrdered(order - 2)).size());
            for (const auto &polyhedronPair : mesh.getPolyhedra())
            {
                std::size_t polyhedronId = polyhedronPair.first;
                const Polyhedron<Polygon3D> &polyhedron = polyhedronPair.second;

                // Get the monomials in 3D up to order k-2, where k order of the VEM
                auto monomials3D = Monomial3D::getMonomialsOrdered(order - 2);

                // Create PolyhedronDof instance for each dof
                for (const auto &monomial : monomials3D)
                {
                    std::shared_ptr<PolyhedronDof> pDof = std::make_shared<PolyhedronDof>(polyhedronId, monomial);
                    dofs.push_back(pDof);
                }
            }
        }

        // You can similarly add other types of dofs (edges, faces, polyhedra) here if needed.
    }

    std::size_t getnumVdofs() const
    {
        return numVdofs;
    }

    std::size_t getnumEdofs() const
    {
        return numEdofs;
    }

    std::size_t getnumFdofs() const
    {
        return numFdofs;
    }

    std::size_t getnumPdofs() const
    {
        return numPdofs;
    }

    // Method to get the corresponding specialized dof to a given id
    template <typename DofType>
    std::shared_ptr<DofType> getDof(std::size_t id) const
    {
        if (id < dofs.size())
        {
            return std::dynamic_pointer_cast<DofType>(dofs[id]);
        }
        else
        {
            return nullptr; // Dof with the given id not found
        }
    }

    // Method to get the corresponding dof to a given id
    std::shared_ptr<VirtualDof> getDof(std::size_t id) const
    {
        if (id < dofs.size())
        {
            return dofs[id];
        }
        else
        {
            return nullptr; // Dof with the given id not found
        }
    }

    // Output stream operator for VirtualDofsCollection
    friend std::ostream &operator<<(std::ostream &os, const VirtualDofsCollection &dofsCollection)
    {
        os << "Virtual degrees of freedom collection - Order: " << order << "\n";
        for (const auto &dof : dofsCollection.dofs)
        {
            dof->operator<<(os);
        }
        return os;
    }
};

// Initialize the static member outside the class definition
unsigned int VirtualDofsCollection::order;

#endif // __VIRTUALDOFS_HPP_