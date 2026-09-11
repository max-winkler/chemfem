#include "fem/FESpace.h"

namespace chemfem{
  namespace fem{

    using chemfem::mesh::Mesh;

    FESpace::FESpace(Mesh& mesh, Element& element, BoundaryIndicator IsDirichlet)
      : refElement(element), mesh(mesh), Dofs(mesh, element, IsDirichlet)
    {
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

    const Element& FESpace::RefElement() const
    {
      return refElement;
    }

    Vector FESpace::IncorporateBC(const Vector& inner) const
    {
      Vector full(Dofs.NrDof());

      for(size_t k=0; k<Dofs.NrDof(); ++k)
	{
	  if(Dofs.IsFree(k))
	    full[k] = inner[Dofs.ReducedIndex(k)];

	  // TODO: Implement also inhomogeneous Dirichlet boundary conditions.
	}

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
	}
      os << " of degree " << space.refElement.Degree() << std::endl;

      os << "Number of DOFs     : " << space.NrDof() << std::endl;
      os << "Number of free DOFs: " << space.NrFreeDof() << std::endl;

      return os;
    }

  }
}
