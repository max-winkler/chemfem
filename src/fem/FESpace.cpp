#include "fem/FESpace.h"

namespace chemfem{
  namespace fem{

    using chemfem::mesh::Cell;
    
    FESpace::FESpace(const Mesh& mesh, Element& element)
      : mesh(mesh), refElement(element)
    {
      if(element.Type() == Lagrange)
	{
	  
	  if(element.Degree() == 1)
	    {
	      DofPerCell = 3;
	      
	      Dof = new size_t[mesh.NrCells()*DofPerCell];

	      int i; std::vector<Cell>::const_iterator cell;
	      for(i=0, cell=mesh.Cells.begin();
		  cell != mesh.Cells.end(); ++cell, ++i)
		{
		  for(int j=0; j<DofPerCell; ++j)
	  	    Dof[i*DofPerCell+j] = cell->LocNode[j]->Index;
		}
	    }
	  else
	    std::cerr << "Error: Lagrange elements of order " << element.Degree() << " not implemented yet\n";
	}
      else
	std::cerr << "Error: Only Lagrange elements are implemented yet.\n";
    }

    size_t FESpace::GetGlobalIndex(size_t cell, size_t index)
    {
      return Dof[DofPerCell*cell + index];
    }

    size_t FESpace::NrDof()
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
