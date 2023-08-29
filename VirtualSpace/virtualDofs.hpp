#ifndef __VIRTUALDOFS_HPP_
#define __VIRTUALDOFS_HPP_

#include "monomial.hpp"
#include "integration.hpp"
#include "mesh.hpp"
#include <memory>
#include <vector>
#include <map>
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
private:
    std::size_t id;
    static std::size_t last_id;

public:
    VirtualDof() : id(last_id++){};
    virtual ~VirtualDof() = default;
    // virtual std::size_t getId() const;
    virtual std::ostream &operator<<(std::ostream &os) const = 0;

    friend std::ostream &operator<<(std::ostream &os, const VertexDof &vDof);
    friend std::ostream &operator<<(std::ostream &os, const EdgeDof &eDof);
    friend std::ostream &operator<<(std::ostream &os, const FaceDof &fDof);
    friend std::ostream &operator<<(std::ostream &os, const PolyhedronDof &pDof);

    virtual std::size_t getId() const
    {
        return id;
    }
};

std::size_t VirtualDof::last_id = 0;

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
        os << "Vertex Dof " << VirtualDof::getId() << " - Vertex ID: " << getId() << ", Coordinates: " << getVertex() << "\n";
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
        os << "Edge Dof " << VirtualDof::getId() << " - Edge ID: " << getId() << ", Gauss-Lobatto Point: " << getGaussLobattoPoint() << ", weight: " << getWeight() << "\n";
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
        os << "Face Dof " << VirtualDof::getId() << " - Face ID: " << getId() << ", Monomial: " << getMonomial() << "\n";
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
        os << "Polyhedron Dof " << VirtualDof::getId() << " - Polyhedron ID: " << getId() << ", Monomial: " << getMonomial() << "\n";
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
            numEdofs = mesh.numEdges() * (order - 1);
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
    }

    static unsigned int getOrder()
    {
        return order;
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

    std::size_t getnumDofs() const
    {
        return dofs.size();
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

    // Print VirtualDofsCollection information
    void print() const
    {
        std::cout<<"VIRTUAL DOFS"<<std::endl;
        std::cout<<"Number of vertex dofs           : "<<this->getnumVdofs()*3<<std::endl;
        std::cout<<"Number of edge dofs             : "<<this->getnumEdofs()*3<<std::endl;
        std::cout<<"Number of face dofs             : "<<this->getnumFdofs()*3<<std::endl;
        std::cout<<"Number of polyhedron dofs       : "<<this->getnumPdofs()*3<<std::endl;
        std::cout<<"Total number of dofs            : "<<this->getnumDofs()*3<<std::endl;
    }
};

// Initialize the static member outside the class definition
unsigned int VirtualDofsCollection::order;

class LocalVirtualDofs
{
private:
    std::vector<std::shared_ptr<VirtualDof>> dofs; // map local dofs to global DOFS
    std::map<std::size_t, std::size_t> V_map;
    std::map<std::size_t, std::vector<std::size_t>> E_map;
    std::map<std::size_t, std::vector<std::size_t>> F_map;
    std::vector<std::size_t> P_vector;

public:
    LocalVirtualDofs(const Polyhedron<Polygon3D> &P, const VirtualDofsCollection &DOFS)
    {
        auto order = DOFS.getOrder();
        std::size_t localDofCounter(0);
        for (std::size_t f = 0; f < P.numPolygons(); f++)
        {
            for (std::size_t e = 0; e < P[f].numEdges(); e++)
            {
                // Insert the vertex into the map and check if it was inserted
                auto inserted = V_map.insert(std::make_pair(P[f][e][0].getId(), localDofCounter));
                // If the vertex was newly inserted, increment local dof counter and add vertex dof to local dofs
                if (inserted.second)
                {
                    ++localDofCounter;
                    dofs.push_back(DOFS.getDof<VertexDof>(P[f][e][0].getId()));
                }
            }
        }
        if (order > 1)
        {
            for (std::size_t f = 0; f < P.numPolygons(); f++)
            {
                for (std::size_t e = 0; e < P[f].numEdges(); e++)
                {
                    //if (P[f][e].getId() > 0)
                    auto EID=P[f].getPositiveEdge(e).getId();
                    if (E_map.find(EID)==E_map.end())
                    {
                        std::vector<std::size_t> Edofs;
                        for (unsigned int i = 0; i < (order - 1); i++)
                        {
                            Edofs.push_back(localDofCounter);
                            auto EDOFId = DOFS.getnumVdofs() + (order - 1) * (EID - 1) + i;
                            // std::cout << EDOFId << " " << std::endl;
                            // std::cout << DOFS.getDof<EdgeDof>(EDOFId)->getGaussLobattoPoint() << std::endl;
                            //  Add edge dof to local dofs
                            dofs.push_back(DOFS.getDof<EdgeDof>(EDOFId));
                            ++localDofCounter;
                        }
                        // Insert the edge dofs into the map
                        E_map.insert(std::make_pair(EID, Edofs));
                        //std::cout<<"inserted a vector of Edofs at position"
                    }
                }
            }
            for (std::size_t f = 0; f < P.numPolygons(); f++)
            {
                auto FId = std::abs(P[f].getId());
                std::vector<std::size_t> Fdofs;
                for (unsigned int i = 0; i < order * (order - 1) / 2; i++)
                {
                    Fdofs.push_back(localDofCounter);
                    auto FDOFId = DOFS.getnumVdofs() + DOFS.getnumEdofs() + (order * (order - 1) / 2) * (FId - 1) + i;
                    // Add face dof to local dofs
                    dofs.push_back(DOFS.getDof<FaceDof>(FDOFId));
                    ++localDofCounter;
                }
                // Insert the face dofs into the map
                F_map.insert(std::make_pair(FId, Fdofs));
            }
            for (unsigned int i = 0; i < (order + 1) * order * (order - 1) / 6; i++)
            {
                P_vector.push_back(localDofCounter);
                auto PDOFId = DOFS.getnumVdofs() + DOFS.getnumEdofs() + DOFS.getnumFdofs() + ((order + 1) * order * (order - 1) / 6) * P.getId() + i;
                // Add polyhedron dof to local dofs
                dofs.push_back(DOFS.getDof<PolyhedronDof>(PDOFId));
                ++localDofCounter;
            }
        }
        /*
        // Print vertices dof
        std::cout << "vertices" << std::endl;
        for (const auto &v : V_map)
        {
            std::cout << v.first << " " << v.second << std::endl;
        }
        for (std::size_t v = 0; v < this->getnumVdofs(); v++)
        {
            std::cout << "Global dof: " << this->getID(v) << " , " << this->getDof<VertexDof>(v)->getVertex() << std::endl;
        }

        // Print edges dof
        std::cout << "edges" << std::endl;
        for (const auto &e : E_map)
        {
            for (const auto &i : e.second)
                std::cout << e.first << " " << i << std::endl;
        }
        for (std::size_t e = this->getnumVdofs(); e < this->getnumVdofs() + this->getnumEdofs(); e++)
        {
            std::cout << "Global dof: " << this->getID(e) << " , " << this->getDof<EdgeDof>(e)->getGaussLobattoPoint() << std::endl;
        }

        // Print faces dof
        std::cout << "faces" << std::endl;
        for (const auto &f : F_map)
        {
            for (const auto &i : f.second)
                std::cout << f.first << " " << i << std::endl;
        }
        for (std::size_t f = this->getnumVdofs() + this->getnumEdofs(); f < this->getnumVdofs() + this->getnumEdofs() + this->getnumFdofs(); f++)
        {
            std::cout << "Global dof: " << this->getID(f) << " , " << this->getDof<FaceDof>(f)->getMonomial() << std::endl;
        }

        // Print polyhedron dof
        std::cout << "polyhedron" << std::endl;
        for (const auto &p : P_vector)
        {
            std::cout << p << std::endl;
        }
        for (std::size_t p = this->getnumVdofs() + this->getnumEdofs() + this->getnumFdofs(); p < this->getnumVdofs() + this->getnumEdofs() + this->getnumFdofs() + this->getnumPdofs(); p++)
        {
            std::cout << "Global dof: " << this->getID(p) << " , " << this->getDof<PolyhedronDof>(p)->getMonomial() << std::endl;
        }
        */
    }

    // Method to get the global id of the corresponding local dof id
    std::size_t getID(std::size_t id) const
    {
        if (id >= dofs.size())
        {
            throw std::out_of_range("Invalid local id."); // Local Dof with the given id does not exist
        }
        return dofs[id]->VirtualDof::getId();
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

    std::size_t VToLocalId(const std::size_t &ID) const { return V_map.at(ID); }
    std::vector<std::size_t> EToLocalId(const std::size_t &ID) const { return E_map.at(ID); }
    std::vector<std::size_t> FToLocalId(const std::size_t &ID) const { return F_map.at(ID); }
    std::size_t PToLocalId(const std::size_t &ID) const { return P_vector[ID]; }

    std::size_t getnumVdofs() const
    {
        return V_map.size();
    }

    std::size_t getnumEdofs() const
    {
        return E_map.size() * (VirtualDofsCollection::getOrder() - 1);
    }

    std::size_t getnumFdofs() const
    {
        return F_map.size() * (VirtualDofsCollection::getOrder() * (VirtualDofsCollection::getOrder() - 1) / 2);
    }

    std::size_t getnumPdofs() const
    {
        return P_vector.size();
    }

    std::size_t getnumDofs() const
    {
        return dofs.size();
    }

    // Output stream operator for LocalVirtualDofs
    friend std::ostream &operator<<(std::ostream &os, const LocalVirtualDofs &dofsCollection)
    {
        os << "Virtual local degrees of freedom collection - Order: " << VirtualDofsCollection::getOrder() << "\n";
        for (const auto &dof : dofsCollection.dofs)
        {
            dof->operator<<(os);
        }
        return os;
    }
};

class LocalVirtualDofsCollection
{
private:
    std::vector<LocalVirtualDofs> Pdofs;

public:
    LocalVirtualDofsCollection(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh, const VirtualDofsCollection &DOFS)
    {
        for (const auto &P : mesh.getPolyhedra())
        {
            Pdofs.emplace_back(LocalVirtualDofs(P.second, DOFS));
            //std::cout<<"inserted polyhedron "<<P.second.getId()<<std::endl;
        }
    }

    const LocalVirtualDofs &getLocalDofs(const std::size_t &Id) const
    {
        //std::cout<<"asked to get dofs of polyherdon "<<Id<<std::endl;
        return Pdofs[Id];
    }

    size_t numLocalDofsCollection() const
    {
        return Pdofs.size();
    }
};

#endif // __VIRTUALDOFS_HPP_