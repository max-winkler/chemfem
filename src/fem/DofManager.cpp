#include "fem/DofManager.h"

#include <iostream>

namespace chemfem{
  namespace fem{

    using chemfem::mesh::Mesh;
    using chemfem::mesh::Node;
    using chemfem::mesh::Cell;
    using chemfem::mesh::Edge;
    using chemfem::mesh::EdgeType;

    DofManager::DofManager(const Mesh& mesh, const Element& element,
			   const BoundaryIndicator& IsDirichlet)
      : nr_local_dof(element.NrDof()), nr_dof(0), nr_free_dof(0)
    {
      if(3*element.DofsPerVertex() + 3*element.DofsPerEdge() + element.DofsInterior()
	 != element.NrDof())
	std::cerr << "Error: The DOFs of the element on the vertices, edges and interior do "
		  << "not add up to its number of local DOFs.\n";

      CreateDofMap(mesh, element);
      MarkDirichletDofs(mesh, element, IsDirichlet);
    }

    void DofManager::CreateDofMap(const Mesh& mesh, const Element& element)
    {
      const size_t nv = element.DofsPerVertex();
      const size_t ne = element.DofsPerEdge();
      const size_t ni = element.DofsInterior();
      const size_t nc = element.NrComponents();

      const size_t edge_offset = mesh.NrNodes()*nv;
      const size_t interior_offset = edge_offset + mesh.NrEdges()*ne;

      nr_dof = interior_offset + mesh.NrCells()*ni;
      DofMap.resize(mesh.NrCells()*nr_local_dof);

      for(size_t c=0; c<mesh.NrCells(); ++c)
	{
	  const Cell& cell = mesh.Cells[c];
	  size_t* local = &DofMap[c*nr_local_dof];

	  for(int v=0; v<3; ++v)
	    for(size_t j=0; j<nv; ++j)
	      *local++ = cell.LocNode[v]*nv + j;

	  for(int k=0; k<3; ++k)
	    {
	      const size_t e = cell.LocEdge[k];

	      // The edge DOFs are numbered from Node0 to Node1 of the edge, a cell that
	      // runs through the edge the other way takes them in reverse order
	      const bool reversed = (mesh.Edges[e].Node0 != cell.LocNode[k]);

	      for(size_t j=0; j<ne; ++j)
		{
		  // Only the scalar DOFs along the edge are mirrored, the components of one
		  // of them keep their order
		  const size_t mirrored = (ne/nc - 1 - j/nc)*nc + j%nc;

		  *local++ = edge_offset + e*ne + (reversed ? mirrored : j);
		}
	    }

	  for(size_t j=0; j<ni; ++j)
	    *local++ = interior_offset + c*ni + j;
	}
    }

    void DofManager::MarkDirichletDofs(const Mesh& mesh, const Element& element,
				       const BoundaryIndicator& IsDirichlet)
    {
      const size_t nv = element.DofsPerVertex();
      const size_t ne = element.DofsPerEdge();
      const size_t edge_offset = mesh.NrNodes()*nv;

      Free.assign(nr_dof, true);
      DirichletEdge.assign(mesh.NrEdges(), false);

      for(size_t e=0; e<mesh.NrEdges(); ++e)
	{
	  const Edge& edge = mesh.Edges[e];

	  if(edge.Type() != EdgeType::BOUNDARY_EDGE)
	    continue;

	  const Node& P0 = mesh.Nodes[edge.Node0];
	  const Node& P1 = mesh.Nodes[edge.Node1];
	  const chemfem::linalg::Coordinate midpoint{0.5*(P0.getX() + P1.getX()),
						     0.5*(P0.getY() + P1.getY())};

	  if(IsDirichlet && !IsDirichlet(midpoint))
	    continue;

	  DirichletEdge[e] = true;

	  for(size_t j=0; j<nv; ++j)
	    {
	      Free[edge.Node0*nv + j] = false;
	      Free[edge.Node1*nv + j] = false;
	    }

	  for(size_t j=0; j<ne; ++j)
	    Free[edge_offset + e*ne + j] = false;
	}

      // Free and Dirichlet DOFs are numbered separately, both in ascending order
      Reduced.resize(nr_dof);

      size_t dirichlet_ctr = 0;
      nr_free_dof = 0;

      for(size_t k=0; k<nr_dof; ++k)
	Reduced[k] = Free[k] ? nr_free_dof++ : dirichlet_ctr++;
    }

    size_t DofManager::NrDof() const
    {
      return nr_dof;
    }

    size_t DofManager::NrFreeDof() const
    {
      return nr_free_dof;
    }

    size_t DofManager::NrLocalDof() const
    {
      return nr_local_dof;
    }

    size_t DofManager::GlobalIndex(size_t cell, size_t local) const
    {
      return DofMap[cell*nr_local_dof + local];
    }

    const size_t* DofManager::LocalDofMap(size_t cell) const
    {
      return &DofMap[cell*nr_local_dof];
    }

    bool DofManager::IsFree(size_t dof) const
    {
      return Free[dof];
    }

    size_t DofManager::ReducedIndex(size_t dof) const
    {
      return Reduced[dof];
    }

    bool DofManager::IsDirichletEdge(size_t edge) const
    {
      return DirichletEdge[edge];
    }

  }
}
