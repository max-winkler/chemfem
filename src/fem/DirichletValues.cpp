#include "fem/DirichletValues.h"

#include <iostream>

#include "mesh/Mesh.h"

using chemfem::linalg::Coordinate;
using chemfem::linalg::Matrix2D;

using chemfem::mesh::Mesh;
using chemfem::mesh::Node;

namespace chemfem{
  namespace fem{

    DirichletValues::DirichletValues(const FESpace& Space)
      : Space(&Space), Values(Space.NrDof(), 0.) {}

    namespace {

      /// Whether values can be prescribed on the space at all
      bool CanPrescribe(const FESpace& Space)
      {
        if(Space.AsVector())
          {
            std::cerr << "Error: The DOFs of this space are fluxes through the edges, not "
                      << "point values, so Dirichlet values cannot be interpolated into "
                      << "it.\n";
            return false;
          }

        if(Space.AllDofsInterior())
          {
            std::cerr << "Error: No DOF of this space touches the boundary, so essential "
                      << "conditions cannot be imposed on it. A DG space takes boundary "
                      << "conditions weakly instead.\n";
            return false;
          }

        return true;
      }
    }

    void DirichletValues::Set(ScalarFunction Value, BoundaryIndicator part)
    {
      if(!CanPrescribe(*Space))
        return;

      const Mesh& mesh = Space->GetMesh();
      const Element& E = Space->RefElement();

      for(size_t c=0; c<mesh.NrCells(); ++c)
	{
	  const Node& x0 = mesh.Nodes[mesh.Cells[c].LocNode[0]];
	  const Coordinate b{x0.getX(), x0.getY()};

	  const Matrix2D Jac = mesh.Jacobian(c);

	  for(size_t k=0; k<Space->NrLocalDof(); ++k)
	    {
	      const size_t dof = Space->GetGlobalIndex(c, k);

	      if(Space->Dofs.IsFree(dof))
		continue;

	      const Coordinate x = b + Jac*E.NodalPoint(k);

	      if(!part || part(x))
		Values[dof] = Value(x);
	    }
	}
    }

    void DirichletValues::Set(VectorFunction Value, BoundaryIndicator part)
    {
      if(!CanPrescribe(*Space))
        return;

      const Mesh& mesh = Space->GetMesh();
      const Element& E = Space->RefElement();

      if(E.NrComponents() != 2)
	{
	  std::cerr << "Error: A vector valued Dirichlet value needs a space with two "
		    << "components.\n";
	  return;
	}

      for(size_t c=0; c<mesh.NrCells(); ++c)
	{
	  const Node& x0 = mesh.Nodes[mesh.Cells[c].LocNode[0]];
	  const Coordinate b{x0.getX(), x0.getY()};

	  const Matrix2D Jac = mesh.Jacobian(c);

	  for(size_t k=0; k<Space->NrLocalDof(); ++k)
	    {
	      const size_t dof = Space->GetGlobalIndex(c, k);

	      if(Space->Dofs.IsFree(dof))
		continue;

	      const Coordinate x = b + Jac*E.NodalPoint(k);

	      if(!part || part(x))
		Values[dof] = Value(x)[E.Component(k)];
	    }
	}
    }

    double DirichletValues::operator[](size_t dof) const
    {
      return Values[dof];
    }

    const FESpace& DirichletValues::GetFESpace() const
    {
      return *Space;
    }

  }
}
