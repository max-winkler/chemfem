#include "fem/FESpace.h"

namespace chemfem{
  namespace fem{

    using chemfem::mesh::Cell;
    using chemfem::mesh::Edge;
    using chemfem::mesh::EdgeType;
    
    FESpace::FESpace(Mesh& mesh, Element& element)
      : refElement(element), mesh(mesh)
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
	      int NrNeighs = edge->Type() == EdgeType::BOUNDARY_EDGE ? 1:2;
	      for(int loc_ind=0; loc_ind < NrNeighs; ++loc_ind)
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
		  if(&(edge->GetNode(0)) == cur_cell.LocNode[edge_ind]
		     && &(edge->GetNode(1)) == cur_cell.LocNode[(edge_ind+1)%3])
		    orientation = true;
		  else if(&(edge->GetNode(1)) == cur_cell.LocNode[edge_ind]
			  && &(edge->GetNode(0)) == cur_cell.LocNode[(edge_ind+1)%3])
		    orientation = false;
		  else
		    std::cerr << "An unexpected error occured. Maybe the mesh data structure "
			      << "is broken.\n";
		  
		  for(int k=0; k<DofPerEdge; ++k)
		    Dof[cell_ind*DofPerCell + 3 + edge_ind*DofPerEdge + k]
		      = NodeDofs + EdgeDofCtr + (orientation ? k : DofPerEdge - k - 1) ;
		}
	      EdgeDofCtr += DofPerEdge;
	    }

	  nr_dof = NodeDofs + EdgeDofCtr;

	  if(element.Degree() >2)
	    {
	      // Numerate interior dofs
	      std::vector<Cell>::const_iterator cell;
	      for(cell = mesh.Cells.begin(); cell != mesh.Cells.end(); ++cell)
		{
		  size_t cell_ind = cell->Index();
		  
		  // TODO: Works only for P3 elements
		  Dof[cell_ind*DofPerCell + 3 + 3*DofPerEdge] = nr_dof++;		  
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
      return nr_dof;
    }

    const Element& FESpace::RefElement() const
    {
      return refElement;
    }
  }
}
