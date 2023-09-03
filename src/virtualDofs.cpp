#include "virtualDofs.hpp"

std::size_t VirtualDof::last_id = 0;

std::size_t VertexDof::getId() const
{
    return vertexId;
}

const Point3D &VertexDof::getVertex() const
{
    return vertex;
}

// Define the stream output operator outside the class definition
std::ostream &operator<<(std::ostream &os, const VertexDof &vDof)
{
    os << "Vertex Dof - Vertex ID: " << vDof.getId() << ", Coordinates: " << vDof.getVertex() << "\n";
    return os;
}

std::size_t EdgeDof::getId() const
{
    return edgeId;
}

const Point3D &EdgeDof::getGaussLobattoPoint() const
{
    return gaussLobattoPoint;
}

const real &EdgeDof::getWeight() const
{
    return weight;
}

// Define the stream output operator outside the class definition
std::ostream &operator<<(std::ostream &os, const EdgeDof &eDof)
{
    os << "Edge Dof - Edge ID: " << eDof.getId() << ", Coordinates: " << eDof.getGaussLobattoPoint() << ", weight: " << eDof.getWeight() << "\n";
    return os;
}

std::size_t FaceDof::getId() const
{
    return faceId;
}

const Monomial2D &FaceDof::getMonomial() const
{
    return monomial;
}

// Define the stream output operator outside the class definition
std::ostream &operator<<(std::ostream &os, const FaceDof &fDof)
{
    os << "Face Dof - Face ID: " << fDof.getId() << ", Monomial: " << fDof.getMonomial() << "\n";
    return os;
}

std::size_t PolyhedronDof::getId() const
{
    return polyhedronId;
}

const Monomial3D &PolyhedronDof::getMonomial() const
{
    return monomial;
}

// Define the stream output operator outside the class definition
std::ostream &operator<<(std::ostream &os, const PolyhedronDof &pDof)
{
    os << "Polyhedron Dof - Polyhedron ID: " << pDof.getId() << ", Monomial: " << pDof.getMonomial() << "\n";
    return os;
}

// Initialize the static member outside the class definition
unsigned int VirtualDofsCollection::order;

VirtualDofsCollection::VirtualDofsCollection(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh, unsigned int order_)
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

unsigned int VirtualDofsCollection::getOrder()
{
    return order;
}

std::size_t VirtualDofsCollection::getnumVdofs() const
{
    return numVdofs;
}

std::size_t VirtualDofsCollection::getnumEdofs() const
{
    return numEdofs;
}

std::size_t VirtualDofsCollection::getnumFdofs() const
{
    return numFdofs;
}

std::size_t VirtualDofsCollection::getnumPdofs() const
{
    return numPdofs;
}

std::size_t VirtualDofsCollection::getnumDofs() const
{
    return dofs.size();
}

std::shared_ptr<VirtualDof> VirtualDofsCollection::getDof(std::size_t id) const
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

std::ostream &operator<<(std::ostream &os, const VirtualDofsCollection &dofsCollection)
{
    os << "Virtual degrees of freedom collection - Order: " << VirtualDofsCollection::order << "\n";
    for (const auto &dof : dofsCollection.dofs)
    {
        dof->operator<<(os);
    }
    return os;
}

void VirtualDofsCollection::print() const
{
    std::cout << "VIRTUAL DOFS" << std::endl;
    std::cout << "Number of vertex dofs           : " << this->getnumVdofs() * 3 << std::endl;
    std::cout << "Number of edge dofs             : " << this->getnumEdofs() * 3 << std::endl;
    std::cout << "Number of face dofs             : " << this->getnumFdofs() * 3 << std::endl;
    std::cout << "Number of polyhedron dofs       : " << this->getnumPdofs() * 3 << std::endl;
    std::cout << "Total number of dofs            : " << this->getnumDofs() * 3 << std::endl;
}

LocalVirtualDofs::LocalVirtualDofs(const Polyhedron<Polygon3D> &P, const VirtualDofsCollection &DOFS)
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
                auto EID = P[f].getPositiveEdge(e).getId();
                if (E_map.find(EID) == E_map.end())
                {
                    std::vector<std::size_t> Edofs;
                    for (unsigned int i = 0; i < (order - 1); i++)
                    {
                        Edofs.push_back(localDofCounter);
                        auto EDOFId = DOFS.getnumVdofs() + (order - 1) * (EID - 1) + i;
                        //  Add edge dof to local dofs
                        dofs.push_back(DOFS.getDof<EdgeDof>(EDOFId));
                        ++localDofCounter;
                    }
                    // Insert the edge dofs into the map
                    E_map.insert(std::make_pair(EID, Edofs));
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
}

std::size_t LocalVirtualDofs::getID(std::size_t id) const
{
    if (id >= dofs.size())
    {
        throw std::out_of_range("Invalid local id."); // Local Dof with the given id does not exist
    }
    return dofs[id]->VirtualDof::getId();
}

std::size_t LocalVirtualDofs::getnumVdofs() const
{
    return V_map.size();
}

std::size_t LocalVirtualDofs::getnumEdofs() const
{
    return E_map.size() * (VirtualDofsCollection::getOrder() - 1);
}

std::size_t LocalVirtualDofs::getnumFdofs() const
{
    return F_map.size() * (VirtualDofsCollection::getOrder() * (VirtualDofsCollection::getOrder() - 1) / 2);
}

std::size_t LocalVirtualDofs::getnumPdofs() const
{
    return P_vector.size();
}

std::size_t LocalVirtualDofs::getnumDofs() const
{
    return dofs.size();
}

std::ostream &operator<<(std::ostream &os, const LocalVirtualDofs &dofsCollection)
{
    os << "Virtual local degrees of freedom collection - Order: " << VirtualDofsCollection::getOrder() << "\n";
    for (const auto &dof : dofsCollection.dofs)
    {
        dof->operator<<(os);
    }
    return os;
}

LocalVirtualDofsCollection::LocalVirtualDofsCollection(const Mesh<Point3D, Edge3D, Polygon3D, Polyhedron<Polygon3D>> &mesh, const VirtualDofsCollection &DOFS)
{
    for (const auto &P : mesh.getPolyhedra())
    {
        Pdofs.emplace_back(LocalVirtualDofs(P.second, DOFS));
        // std::cout<<"inserted polyhedron "<<P.second.getId()<<std::endl;
    }
}

const LocalVirtualDofs &LocalVirtualDofsCollection::getLocalDofs(const std::size_t &Id) const
{
    // std::cout<<"asked to get dofs of polyherdon "<<Id<<std::endl;
    return Pdofs[Id];
}

size_t LocalVirtualDofsCollection::numLocalDofsCollection() const
{
    return Pdofs.size();
}