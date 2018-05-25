#include "fem/FESpace.h"

namespace chemfem{
  namespace fem{

    using chemfem::mesh::Cell;
    
    FESpace::FESpace(const Mesh& mesh, FEType type, int degree)
      : mesh(mesh), type(type), Degree(degree)
    {
      if(type == Lagrange)
	{
	  
	  if(degree == 1)
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
	    std::cerr << "Error: Lagrange elements of order " << degree << " not implemented yet\n";
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
      if(type == Lagrange && Degree == 1)
	return mesh.NrNodes();
      else
	return -1;
    }
    
  }
}
