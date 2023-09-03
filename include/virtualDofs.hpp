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
    /**
     * @brief Constructor
     * 
     */
    VirtualDof() : id(last_id++){};

    /**
     * @brief Destructor
     * 
     */
    virtual ~VirtualDof() = default;

    /**
     * @brief Pure virtual output stream operator
     * 
     * @param os 
     * @return std::ostream& 
     */
    virtual std::ostream &operator<<(std::ostream &os) const = 0;

    /**
     * @brief Output stream operators
     * 
     * @param os 
     * @param vDof 
     * @return std::ostream& 
     */
    friend std::ostream &operator<<(std::ostream &os, const VertexDof &vDof);
    friend std::ostream &operator<<(std::ostream &os, const EdgeDof &eDof);
    friend std::ostream &operator<<(std::ostream &os, const FaceDof &fDof);
    friend std::ostream &operator<<(std::ostream &os, const PolyhedronDof &pDof);

    /**
     * @brief Getter for the VirtualDof Id
     * 
     * @return std::size_t 
     */
    virtual std::size_t getId() const
    {
        return id;
    }
};

// Derived class for VertexDof
class VertexDof : public VirtualDof
{
private:
    std::size_t vertexId;
    Point3D vertex;

public:
    /**
     * @brief Constructor
     * 
     * @param id 
     * @param vertex_ 
     */
    VertexDof(std::size_t id, const Point3D &vertex_) : vertexId(id), vertex(vertex_) {}

    /**
     * @brief Getter fot the id of the vertex
     * 
     * @return std::size_t 
     */
    std::size_t getId() const override;

    /**
     * @brief Getter for the vertex
     * 
     * @return const Point3D& 
     */
    const Point3D &getVertex() const;

    /**
     * @brief Output stream operator
     * 
     * @param os 
     * @return std::ostream& 
     */
    std::ostream &operator<<(std::ostream &os) const override
    {
        os << "Vertex Dof " << VirtualDof::getId() << " - Vertex ID: " << getId() << ", Coordinates: " << getVertex() << "\n";
        return os;
    }
};

// Derived class for EdgeDof
class EdgeDof : public VirtualDof
{
private:
    std::size_t edgeId;
    Point3D gaussLobattoPoint;
    real weight;

public:
    /**
     * @brief Constructor
     * 
     * @param id 
     * @param gaussLobattoPoint_ 
     * @param weight_ 
     */
    EdgeDof(std::size_t id, const Point3D &gaussLobattoPoint_, const real &weight_)
        : edgeId(id), gaussLobattoPoint(gaussLobattoPoint_), weight(weight_) {}

    /**
     * @brief Getter for the id of the Edge DOF
     * 
     * @return std::size_t 
     */
    std::size_t getId() const override;

    /**
     * @brief Getter for the Gauss-Lobatto point
     * 
     * @return const Point3D& 
     */
    const Point3D &getGaussLobattoPoint() const;

    /**
     * @brief Getter for the weight of the Gauss-Lobatto point
     * 
     * @return const real& 
     */
    const real &getWeight() const;

    /**
     * @brief Output stream operator
     * 
     * @param os 
     * @return std::ostream& 
     */
    std::ostream &operator<<(std::ostream &os) const override
    {
        os << "Edge Dof " << VirtualDof::getId() << " - Edge ID: " << getId() << ", Gauss-Lobatto Point: " << getGaussLobattoPoint() << ", weight: " << getWeight() << "\n";
        return os;
    }
};

// Derived class for FaceDof
class FaceDof : public VirtualDof
{
private:
    std::size_t faceId;
    Monomial2D monomial;

public:
    /**
     * @brief Constructor
     * 
     * @param id 
     * @param monomial_ 
     */
    FaceDof(std::size_t id, const Monomial2D &monomial_) : faceId(id), monomial(monomial_) {}

    /**
     * @brief Getter for the id of the face
     * 
     * @return std::size_t 
     */
    std::size_t getId() const override;

    /**
     * @brief Getter for the monomial of the Face DOF
     * 
     * @return const Monomial2D& 
     */
    const Monomial2D &getMonomial() const;

    /**
     * @brief Output stream operator
     * 
     * @param os 
     * @return std::ostream& 
     */
    std::ostream &operator<<(std::ostream &os) const override
    {
        os << "Face Dof " << VirtualDof::getId() << " - Face ID: " << getId() << ", Monomial: " << getMonomial() << "\n";
        return os;
    }
};

// Derived class for PolyhedronDof
class PolyhedronDof : public VirtualDof
{
private:
    std::size_t polyhedronId;
    Monomial3D monomial;

public:
    /**
     * @brief Constructor
     * 
     * @param id 
     * @param monomial_ 
     */
    PolyhedronDof(std::size_t id, const Monomial3D &monomial_) : polyhedronId(id), monomial(monomial_) {}

    /**
     * @brief Getter for the id of the polyhedron
     * 
     * @return std::size_t 
     */
    std::size_t getId() const override;

    /**
     * @brief Getter for the monomial of the polyhedron DOF
     * 
     * @return const Monomial3D& 
     */
    const Monomial3D &getMonomial() const;

    /**
     * @brief Output stream operator
     * 
     * @param os 
     * @return std::ostream& 
     */
    std::ostream &operator<<(std::ostream &os) const override
    {
        os << "Polyhedron Dof " << VirtualDof::getId() << " - Polyhedron ID: " << getId() << ", Monomial: " << getMonomial() << "\n";
        return os;
    }
};

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
    /**
     * @brief Virtual default destructor
     * 
     */
    virtual ~VirtualDofsCollection() = default;

    /**
     * @brief Constructor to create VirtualDofs from a Mesh and order
     * 
     * @param mesh 
     * @param order_ 
     */
    VirtualDofsCollection(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh, unsigned int order_);

    /**
     * @brief Get the order of the virtual dofs collection
     * 
     * @return unsigned int 
     */
    static unsigned int getOrder();

    /**
     * @brief Get the number of vertex-type dofs
     * 
     * @return std::size_t 
     */
    std::size_t getnumVdofs() const;

    /**
     * @brief Get the number of edge-type dofs
     * 
     * @return std::size_t 
     */
    std::size_t getnumEdofs() const;

    /**
     * @brief Get the number of face-type dofs
     * 
     * @return std::size_t 
     */
    std::size_t getnumFdofs() const;

    /**
     * @brief Get the number of polyhedron-type dofs
     * 
     * @return std::size_t 
     */
    std::size_t getnumPdofs() const;

    /**
     * @brief Get the total number of dofs
     * 
     * @return std::size_t 
     */
    std::size_t getnumDofs() const;

    /**
     * @brief Method to get the corresponding specialized dof to a given id
     * 
     * @tparam DofType 
     * @param id 
     * @return std::shared_ptr<DofType> 
     */
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

    /**
     * @brief Method to get the corresponding dof to a given id
     * 
     * @param id 
     * @return std::shared_ptr<VirtualDof> 
     */
    std::shared_ptr<VirtualDof> getDof(std::size_t id) const;

    /**
     * @brief Output stream operator for VirtualDofsCollection
     * 
     * @param os 
     * @param dofsCollection 
     * @return std::ostream& 
     */
    friend std::ostream &operator<<(std::ostream &os, const VirtualDofsCollection &dofsCollection);

    /**
     * @brief Print VirtualDofsCollection information
     * 
     */
    void print() const;
};

class LocalVirtualDofs
{
private:
    std::vector<std::shared_ptr<VirtualDof>> dofs; // map local dofs to global DOFS
    std::map<std::size_t, std::size_t> V_map;
    std::map<std::size_t, std::vector<std::size_t>> E_map;
    std::map<std::size_t, std::vector<std::size_t>> F_map;
    std::vector<std::size_t> P_vector;

public:
    /**
     * @brief Constructor
     * 
     * @param P 
     * @param DOFS 
     */
    LocalVirtualDofs(const Polyhedron<Polygon3D> &P, const VirtualDofsCollection &DOFS);

    /**
     * @brief Method to get the global id of the corresponding local dof id
     * 
     * @param id 
     * @return std::size_t 
     */
    std::size_t getID(std::size_t id) const;

    /**
     * @brief Method to get the corresponding specialized dof to a given id
     * 
     * @tparam DofType 
     * @param id 
     * @return std::shared_ptr<DofType> 
     */
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

    /**
     * @brief Get the corresponding local dof for vertex-type global dof
     * 
     * @param ID 
     * @return std::size_t 
     */
    std::size_t VToLocalId(const std::size_t &ID) const { return V_map.at(ID); }

    /**
     * @brief Get the corresponding local dof for edge-type global dof
     * 
     * @param ID 
     * @return std::vector<std::size_t> 
     */
    std::vector<std::size_t> EToLocalId(const std::size_t &ID) const { return E_map.at(ID); }

    /**
     * @brief Get the corresponding local dof for face-type global dof
     * 
     * @param ID 
     * @return std::vector<std::size_t> 
     */
    std::vector<std::size_t> FToLocalId(const std::size_t &ID) const { return F_map.at(ID); }

    /**
     * @brief Get the corresponding local dof for polyhedron-type global dof
     * 
     * @param ID 
     * @return std::size_t 
     */
    std::size_t PToLocalId(const std::size_t &ID) const { return P_vector[ID]; }

    /**
     * @brief Get the number of vertex-type dofs
     * 
     * @return std::size_t 
     */
    std::size_t getnumVdofs() const;

    /**
     * @brief Get the number of edge-type dofs
     * 
     * @return std::size_t 
     */
    std::size_t getnumEdofs() const;

    /**
     * @brief Get the number of face-type dofs
     * 
     * @return std::size_t 
     */
    std::size_t getnumFdofs() const;

    /**
     * @brief Get the number of polyhedron-type dofs
     * 
     * @return std::size_t 
     */
    std::size_t getnumPdofs() const;

    /**
     * @brief Get the total number of local dofs
     * 
     * @return std::size_t 
     */
    std::size_t getnumDofs() const;

    /**
     * @brief Output stream operator for LocalVirtualDofs
     * 
     * @param os 
     * @param dofsCollection 
     * @return std::ostream& 
     */
    friend std::ostream &operator<<(std::ostream &os, const LocalVirtualDofs &dofsCollection);
};

class LocalVirtualDofsCollection
{
private:
    std::vector<LocalVirtualDofs> Pdofs;

public:
    /**
     * @brief Constructor
     * 
     * @param mesh 
     * @param DOFS 
     */
    LocalVirtualDofsCollection(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh, const VirtualDofsCollection &DOFS);

    /**
     * @brief Get the LocalVirtualDofs of the element Id
     * 
     * @param Id 
     * @return const LocalVirtualDofs& 
     */
    const LocalVirtualDofs &getLocalDofs(const std::size_t &Id) const;

    /**
     * @brief Get the number of LocalVirtualDofs in the collection
     * 
     * @return size_t 
     */
    size_t numLocalDofsCollection() const;
};

#endif // __VIRTUALDOFS_HPP_