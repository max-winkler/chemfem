#include "fem/FESpace.h"

namespace chemfem{
  namespace fem{

    using chemfem::mesh::Cell;
    
    FESpace::FESpace(const Mesh& mesh, Element& element)
      : mesh(mesh), refElement(element)
    {
      DofPerCell = element.NrDof();
      Dof = new size_t[mesh.NrCells()*DofPerCell];
      
      if(element.Type() == Lagrange)
	{
	  // Dofs in the vertices
	  int i; std::vector<Cell>::const_iterator cell;
	  for(i=0, cell=mesh.Cells.begin();
	      cell != mesh.Cells.end(); ++cell, ++i)
	    {
	      for(int j=0; j<3; ++j)
		Dof[i*DofPerCell+j] = cell->LocNode[j]->Index;
	    }

	  size_t NodeDofs = mesh.NrNodes();
	  
	  // Dofs at the edges
	  
	}      
      else
	std::cerr << "Error: Only Lagrange elements are implemented yet.\n";
    }

    size_t FESpace::GetGlobalIndex(size_t cell, size_t index) const
    {
      return Dof[DofPerCell*cell + index];
    }

    size_t FESpace::NrDof() const
    {
      if(refElement.Type() == Lagrange && refElement.Degree() == 1)
	return mesh.NrNodes();
      else
	return -1;
    }

    const Element& FESpace::RefElement() const
    {
      return refElement;
    }
  }
}
