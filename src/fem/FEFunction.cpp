#include "fem/FEFunction.h"

namespace chemfem{
  namespace fem{

    using chemfem::mesh::Cell;
    using chemfem::mesh::Node;
    using chemfem::mesh::Edge;
    
    FEFunction::FEFunction(const FEFunction& other)
      : Space(other.Space), Data(other.Data) {}

    FEFunction::FEFunction(const FESpace& Space)
      : Space(&Space), Data(Space.NrDof()) {}

    const FESpace& FEFunction::GetFESpace() const
    {
      return *Space;
    }

    const double& FEFunction::operator[](size_t k) const
    {
      return Data[k];
    }

    void FEFunction::CreateFunction(const Vector& FreeDof)
    {
      Data = Space->IncorporateBC(FreeDof);
    }

    void FEFunction::SetCoefficients(const Vector& Data)
    {
      this->Data = Data;
    }

    FEFunction FESpace::Interpolate(double (*u)(double, double))
    {
      // TODO: This functions seems to be wrong for P2 elements.
      // Check this!
      FEFunction Function(*this);
      Vector Vec(nr_dof);

      double x, y;
      
      if(refElement.Type() == Lagrange)
	{
	  // Interpolate nodes
	  std::vector<Cell>::const_iterator it_cell;
	  for(it_cell = mesh.Cells.begin(); it_cell != mesh.Cells.end(); ++it_cell)
	    {
	      const Cell& cell = *it_cell;
	      size_t cell_ind = cell.Index();

	      const size_t* LocalDof = GetLocalDofMap(cell_ind);
	      
	      // Interpolate in nodes of the mesh
	      for(int k=0; k<3; ++k)
		{
		  Node& node = mesh.Nodes[cell.LocNode[k]];
		  x = node.getX();
		  y = node.getY();

		  // TODO: Each cell sharing this note has written this variable already.
		  // This may lead to performance lost. Can we improve this?
		  Vec[LocalDof[k]] = u(x,y);
		}

	      if(refElement.Degree() < 2) continue;
	      
	      // Interpolate in nodes at edges
	      for(int l=0; l<3; ++l)
		{
		  const Edge& edge = mesh.Edges[cell.LocEdge[l]];

		  // Coordinates of endpoints
		  double x0 = mesh.Nodes[edge.Node0].getX(), y0 = mesh.Nodes[edge.Node0].getY();
		  double x1 = mesh.Nodes[edge.Node1].getX(), y1 = mesh.Nodes[edge.Node1].getY();
		  
		  // Determine orientation of the edge
		  bool orientation = true;
		  if(edge.Node0 == cell.LocNode[l] && edge.Node1 == cell.LocNode[(l+1)%3])
		    orientation = true;
		  else if(edge.Node1 == cell.LocNode[l] && edge.Node0 == cell.LocNode[(l+1)%3])
		    orientation = false;
		  else
		    std::cerr << "An unexpected error occured. Maybe the mesh data structure "
			      << "is broken.\n";
		  
		  for(int k=0; k<DofPerEdge; ++k)
		    {
		      // Coordinate of edge dof
		      if(orientation)
			{
			  x = x0 + (x1-x0)*((double)k/DofPerEdge);
			  y = y0 + (y1-y0)*((double)k/DofPerEdge);
			}
		      else
			{
			  x = x1 + (x0-x1)*((double)k/DofPerEdge);
			  y = y1 + (y0-y1)*((double)k/DofPerEdge);
			}
		      
		      Vec[DofMap[cell_ind*DofPerCell + 3 + l*DofPerEdge + k]]
			= u(x,y);		      
		    } // loop over edges

		  if(refElement.Degree() < 3) continue;

		  double X[3], Y[3];
		  for(int k=0; k<3; ++k)
		    {
		      X[k] = mesh.Nodes[cell.LocNode[k]].getX();
		      Y[k] = mesh.Nodes[cell.LocNode[k]].getY();
		    }
		  
		  // Interpolate interior nodes		  
		  for(int k=0; k<IntDofPerCell; ++k)
		    {
		      // TODO: Is there a more elegant way to determine interior DOF's?
		      // Now we are restricted to P3 and P4 elements.
		      switch(refElement.Degree())
			{
			case 3:
			  x = (X[0]+X[1]+X[2])/3;
			  y = (Y[0]+Y[1]+Y[2])/3;
			  break;
			case 4:
			  x = (X[0]+X[1]+X[2]+X[k])/4;
			  y = (Y[0]+Y[1]+Y[2]+Y[k])/4;
			  break;
			default:
			  std::cerr << "Interpolation of P" << refElement.Degree()
				    << " not implemented yet.\n"; 
			  break;
			}
		      Vec[DofMap[cell_ind*DofPerCell + 3 + 3*DofPerEdge + k]]
			= u(x,y);
		    }
		  
		} // loop over cells
	    }
	  
	  Function.SetCoefficients(Vec);
	}
      else
	std::cerr << "Only Lagrange elements implemented. Cannot interpolate yet.\n";

      return Function;
    }

    void FEFunction::WriteVtk(const std::string& filename) const
    {
      Space->GetMesh().WriteVtk(filename, Data);
    }
  }
}
