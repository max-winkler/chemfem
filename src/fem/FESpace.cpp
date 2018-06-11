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
	  size_t EdgeDofCtr = 0;
	  // Dofs at the edges	  
	  std::set<Edge>::iterator edge;
	  for(edge = mesh.Edges.begin(); edge != mesh.Edges.end(); ++edge)
	    {
	      for(int loc_ind=0; loc_ind < 2; ++loc_ind)
		{
		  Cell& cur_cell = edge->GetNeighbor(loc_ind);
		  int cell_ind = cur_cell.Index();

		  int edge_ind;
		  for(edge_ind=0; edge_ind<3; ++edge_ind)
		    if(cur_cell.LocEdge[edge_ind] == &(*edge)) break;
		  if(edge_ind == 3)
		    {
		      std::cerr << "Edge not found but it should belong to the element. "
				<< "Maybe the mesh format is corrupt.\n";
		      return;
		    }

		  bool orientation;
		  if(edge->Node0 == cur_cell.LocNode[edge_ind]
		     && edge->Node1 == cur_cell.LocNode[(edge_ind+1)%3])
		    orientation = true;
		  else if(edge->Node1 == cur_cell.LocNode[edge_ind]
			  && edge->Node0 == cur_cell.LocNode[(edge_ind+1)%3])
		    orientation = false;
		  else
		    std::cerr << "An unexpected error occured. Maybe the mesh format is corrupt.\n";
		  
		  for(int k=0; k<DofPerEdge; ++k)
		    Dof[cell_ind*DofPerCell + 3 + k] = NodeDofs + EdgeDofCtr + k;
		  // TODO: DOF depends on local edge index and orientation.
		}
	      EdgeDofCtr += DofPerEdge;
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
