#include "fem/FESpace.h"

namespace chemfem{
  namespace fem{

    using chemfem::mesh::Cell;
    using chemfem::mesh::Edge;
    
    FESpace::FESpace(Mesh& mesh, Element& element)
      : mesh(mesh), refElement(element)
    {
      DofPerCell = element.NrDof();
      DofPerEdge = element.Degree() - 1;
	
      Dof = new size_t[mesh.NrCells()*DofPerCell];
      
      if(element.Type() == Lagrange)
	{
	  // Dofs in the vertices
	  size_t i; std::vector<Cell>::iterator cell;
	  for(i=0, cell = mesh.Cells.begin();
	      cell != mesh.Cells.end(); ++cell, ++i)
	    {
	      for(int j=0; j<3; ++j)
		Dof[i*DofPerCell+j] = cell->LocNode[j]->Index();

	      cell->SetIndex(i);
	    }

	  size_t NodeDofs = mesh.NrNodes();
	  
	  // Dofs at the edges	  
	  std::set<Edge>::iterator edge;
	  for(edge = mesh.Edges.begin(); edge != mesh.Edges.end(); ++edge)
	    {
	      for(int loc_ind=0; loc_ind < 2; ++loc_ind)
		{
		  Cell& cur_cell = edge->GetNeighbor(loc_ind);
		  int cell_ind = cur_cell.Index();

		  for(int k=0; k<DofPerEdge; ++k)
		    Dof[cell_ind*DofPerCell + 3 + k] = 10000; // TODO: Continue here.
		}
	    }	  
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
