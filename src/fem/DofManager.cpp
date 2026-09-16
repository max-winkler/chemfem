#include "fem/DofManager.h"

#include <cmath>
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
      Reversed.assign(3*mesh.NrCells(), false);

      // An edge DOF of a vector element is a flux through a globally fixed normal
      const bool piola = (dynamic_cast<const VectorElement*>(&element) != nullptr);

      if(piola)
	Sign.assign(mesh.NrCells()*nr_local_dof, 1.);

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

	      Reversed[3*c + k] = reversed;

	      for(size_t j=0; j<ne; ++j)
		{
		  // Only the scalar DOFs along the edge are mirrored, the components of one
		  // of them keep their order
		  const size_t mirrored = (ne/nc - 1 - j/nc)*nc + j%nc;

		  *local++ = edge_offset + e*ne + (reversed ? mirrored : j);

		  if(piola && reversed)
		    Sign[c*nr_local_dof + 3*nv + k*ne + j] = -1.;
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
      Frames.assign(mesh.NrNodes(), NodeFrame{Vector2D(1., 0.), Vector2D(0., 1.)});

      // Which DOF of a vertex is its value and which are its derivatives. The derivatives
      // are numbered like the frame directions they are taken along.
      int value_dof = -1, deriv_dof[2] = {-1, -1};
      int nr_deriv = 0;

      for(size_t j=0; j<nv; ++j)
	{
	  const DofDescriptor d = element.Dof(int(j));

	  if(d.type == DofType::PointValue && value_dof < 0)
	    value_dof = int(j);
	  else if(d.type == DofType::EdgeDirectionalDerivative && nr_deriv < 2)
	    deriv_dof[nr_deriv++] = int(j);
	}

      // The restriction of such an element to an edge is determined by the value and the
      // derivative along that edge at its two endpoints, so u = 0 on a straight side leaves
      // the normal derivative free. Without derivative DOFs every DOF of the vertex is fixed.
      const bool turn = (nr_deriv == 2 && value_dof >= 0);

      std::vector<Vector2D> tangent(mesh.NrNodes(), Vector2D(0., 0.));
      std::vector<int> rank(mesh.NrNodes(), 0);

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

	  for(size_t j=0; j<ne; ++j)
	    Free[edge_offset + e*ne + j] = false;

	  if(!turn)
	    {
	      for(size_t j=0; j<nv; ++j)
		{
		  Free[edge.Node0*nv + j] = false;
		  Free[edge.Node1*nv + j] = false;
		}

	      continue;
	    }

	  Vector2D t(P1.getX() - P0.getX(), P1.getY() - P0.getY());
	  const double len = t.Norm();

	  if(len > 0.)
	    t *= 1./len;

	  const size_t ends[2] = {edge.Node0, edge.Node1};

	  for(int i=0; i<2; ++i)
	    {
	      const size_t n = ends[i];

	      // Two sides that are not parallel fix the whole gradient. The tolerance is
	      // tight on purpose: a shallow corner counts as a corner, which is the
	      // condition the exact solution on the polygon satisfies there.
	      if(rank[n] == 0)
		{
		  tangent[n] = t;
		  rank[n] = 1;
		}
	      else if(rank[n] == 1
		      && std::fabs(tangent[n].x*t.y - tangent[n].y*t.x) > 1.e-12)
		rank[n] = 2;
	    }
	}

      if(turn)
	for(size_t n=0; n<mesh.NrNodes(); ++n)
	  {
	    if(rank[n] == 0)
	      continue;

	    Free[n*nv + value_dof] = false;
	    Free[n*nv + deriv_dof[0]] = false;

	    if(rank[n] == 2)
	      {
		Free[n*nv + deriv_dof[1]] = false;
		continue;
	      }

	    Frames[n] = NodeFrame{tangent[n], Vector2D(-tangent[n].y, tangent[n].x)};
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

    double DofManager::LocalSign(size_t cell, size_t local) const
    {
      return Sign.empty() ? 1. : Sign[cell*nr_local_dof + local];
    }

    bool DofManager::IsEdgeReversed(size_t cell, int local_edge) const
    {
      return Reversed[3*cell + local_edge];
    }

    const NodeFrame& DofManager::Frame(size_t node) const
    {
      return Frames[node];
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
