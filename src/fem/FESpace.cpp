#include "fem/FESpace.h"
#include "fem/DirichletValues.h"

namespace chemfem{
  namespace fem{

    using chemfem::mesh::Mesh;

    FESpace::FESpace(Mesh& mesh, Element& element, BoundaryIndicator IsDirichlet)
      : refElement(element), scalar(dynamic_cast<const ScalarElement*>(&element)),
	vector(dynamic_cast<const VectorElement*>(&element)),
	mesh(mesh), Dofs(mesh, element, IsDirichlet)
    {
      if(!scalar && !vector)
	std::cerr << "Error: The element is neither a ScalarElement nor a VectorElement.\n";

      // Test if vertices are numbered correctly
      for(size_t i=0; i<mesh.Nodes.size(); ++i)
	{
	  if(mesh.Nodes[i].Index() != i)
	    std::cerr << "Nodes are not numbered correctly\n";
	}
    }

    size_t FESpace::GetGlobalIndex(size_t cell, size_t index) const
    {
      return Dofs.GlobalIndex(cell, index);
    }

    size_t FESpace::NrDof() const
    {
      return Dofs.NrDof();
    }

    size_t FESpace::NrFreeDof() const
    {
      return Dofs.NrFreeDof();
    }

    size_t FESpace::NrLocalDof() const
    {
      return Dofs.NrLocalDof();
    }

    size_t FESpace::NrComponents() const
    {
      return refElement.NrComponents();
    }

    const Element& FESpace::RefElement() const
    {
      return refElement;
    }

    Vector FESpace::IncorporateBC(const Vector& inner) const
    {
      Vector full(Dofs.NrDof());

      for(size_t k=0; k<Dofs.NrDof(); ++k)
	if(Dofs.IsFree(k))
	  full[k] = inner[Dofs.ReducedIndex(k)];

      return full;
    }

    Vector FESpace::IncorporateBC(const Vector& inner, const DirichletValues& g) const
    {
      Vector full(Dofs.NrDof());

      for(size_t k=0; k<Dofs.NrDof(); ++k)
	full[k] = Dofs.IsFree(k) ? inner[Dofs.ReducedIndex(k)] : g[k];

      return full;
    }

    const Mesh& FESpace::GetMesh() const
    {
      return mesh;
    }

    const size_t* FESpace::GetLocalDofMap(size_t k) const
    {
      return Dofs.LocalDofMap(k);
    }

    double FESpace::LocalSign(size_t cell, size_t local) const
    {
      return Dofs.LocalSign(cell, local);
    }

    bool FESpace::IsEdgeReversed(size_t cell, int local_edge) const
    {
      return Dofs.IsEdgeReversed(cell, local_edge);
    }

    const ScalarElement* FESpace::AsScalar() const
    {
      return scalar;
    }

    const VectorElement* FESpace::AsVector() const
    {
      return vector;
    }

    bool FESpace::IsVectorValued() const
    {
      return vector != nullptr || NrComponents() > 1;
    }

    bool FESpace::AllDofsInterior() const
    {
      return refElement.DofsPerVertex() == 0 && refElement.DofsPerEdge() == 0;
    }

    std::ostream& operator<<(std::ostream& os, const FESpace& space)
    {
      os << "FE Space of type ";
      switch(space.refElement.Type())
	{
	case Lagrange:
	  os << "Lagrange";
	  break;
	case CrouzeixRaviart:
	  os << "Crouzeix-Raviart";
	  break;
	case DG:
	  os << "discontinuous Galerkin";
	  break;
	case RaviartThomas:
	  os << "Raviart-Thomas";
	  break;
	}
      os << " of degree " << space.refElement.Degree() << std::endl;

      os << "Number of DOFs     : " << space.NrDof() << std::endl;
      os << "Number of free DOFs: " << space.NrFreeDof() << std::endl;

      return os;
    }

  }
}
