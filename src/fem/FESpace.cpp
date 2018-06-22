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
      IntDofPerCell = (element.Degree()-1) * (element.Degree()-2) / 2;
      
      CreateDofMap();
      CreateDirichletBcMap();
    }

    void FESpace::CreateDofMap()
    {
      DofMap = new size_t[mesh.NrCells()*DofPerCell];
      
      if(refElement.Type() == Lagrange)
	{
	  // Dofs in the vertices
	  size_t i; std::vector<Cell>::iterator it_cell;
	  for(i=0, it_cell = mesh.Cells.begin();
	      it_cell != mesh.Cells.end(); ++it_cell, ++i)
	    {
	      Cell& cell = *it_cell;
	      
	      for(int j=0; j<3; ++j)
		DofMap[i*DofPerCell+j] = cell.LocNode[j]->Index();

	      cell.SetIndex(i);
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
		  int edge_ind = cur_cell.EdgeIndex(*edge);
		  
		  if(edge_ind == -1)
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
		    DofMap[cell_ind*DofPerCell + 3 + edge_ind*DofPerEdge + k]
		      = NodeDofs + EdgeDofCtr + (orientation ? k : DofPerEdge - k - 1) ;
		}
	      EdgeDofCtr += DofPerEdge;
	    }

	  nr_dof = NodeDofs + EdgeDofCtr;

	  if(refElement.Degree() >2)
	    {
	      // Numerate interior dofs
	      std::vector<Cell>::const_iterator cell;
	      for(cell = mesh.Cells.begin(); cell != mesh.Cells.end(); ++cell)
		{
		  size_t cell_ind = cell->Index();
		  
		  for(int k=0; k<IntDofPerCell; ++k)
		    DofMap[cell_ind*DofPerCell + 3 + 3*DofPerEdge + k] = nr_dof++;		  
		}
	    }
	}
      else
	std::cerr << "Error: Only Lagrange elements are implemented yet.\n";
    }

    void FESpace::CreateDirichletBcMap()
    {      
      std::vector<const Edge*>::const_iterator it_edge;
      for(it_edge = mesh.BdEdges.begin(); it_edge != mesh.BdEdges.end(); ++it_edge)
	{
	  const Edge& edge = *(*it_edge);
	  Cell& cell = edge.GetNeighbor(0);
	  int edge_ind = cell.EdgeIndex(edge);
	  
	  for(int k=0; k<2; ++k)
	    DirichletNodes.insert(edge.GetNode(k).Index());
	  for(int k=0; k<DofPerEdge; ++k)
	    DirichletNodes.insert(DofMap[cell.Index()*DofPerCell + 3 + edge_ind*DofPerEdge + k]);
	} // loop over boundary edges

      DofType = new bool[nr_dof];
      DofIndex = new size_t[nr_dof];
      
      // Initialize DOF type and DOF index arrays
      std::set<size_t>::iterator it_node;
      size_t dirichlet_dof_ctr, free_dof_ctr = 0, last_dirichlet_dof = 0;
      
      for(it_node = DirichletNodes.begin(), dirichlet_dof_ctr=0; it_node != DirichletNodes.end();
	  ++it_node, ++dirichlet_dof_ctr)
	{
	  const size_t index = *(it_node);

	  // Numerate free nodes between last and current Dirichlet Dof
	  for(size_t k = last_dirichlet_dof+1; k<index; ++k)
	    {
	      DofType[k] = true;
	      DofIndex[k] = free_dof_ctr++;
	    }
	  
	  // Numerate Dirichlet Dof
	  DofType[index] = false;
	  DofIndex[index] = dirichlet_dof_ctr;

	  last_dirichlet_dof = index;	  
	}
      for(size_t k = last_dirichlet_dof+1; k<nr_dof; ++k)
	DofIndex[k] = free_dof_ctr++;

      nr_free_dof = free_dof_ctr;
    }
    
    size_t FESpace::GetGlobalIndex(size_t cell, size_t index) const
    {
      return DofMap[DofPerCell*cell + index];
    }

    size_t FESpace::NrDof() const
    {
      return nr_dof;
    }

    size_t FESpace::NrFreeDof() const
    {
      return nr_free_dof;
    }

    const Element& FESpace::RefElement() const
    {
      return refElement;
    }

    Vector FESpace::IncorporateBC(const Vector& inner) const
    {
      Vector full(nr_dof);

      for(size_t k=0; k<nr_dof; ++k)
	{
	  if(DofType[k])
	    full[k] = inner[DofIndex[k]];

	  // TODO: Implement also inhomogeneous Dirichlet boundary conditions.
	}

      return full;
    }
    
  }
}
