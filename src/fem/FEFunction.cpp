#include "fem/FEFunction.h"

namespace chemfem{
  namespace fem{

    using chemfem::mesh::Node;

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
      CacheValid = false;
    }

    void FEFunction::SetCoefficients(const Vector& Data)
    {
      this->Data = Data;
      CacheValid = false;
    }

    PointValues FEFunction::Evaluate(const QuadPoint& p) const
    {
      if(CacheValid && p.cell == CachedPoint.cell
	 && p.xi == CachedPoint.xi && p.eta == CachedPoint.eta)
	return CachedValues;

      const Element& E = Space->RefElement();

      PointValues ref;
      ref.value = 0.;

      for(size_t k=0; k<Space->NrLocalDof(); ++k)
	{
	  const double coeff = Data[Space->GetGlobalIndex(p.cell, k)];

	  ref.value += coeff * E.Value(k, p.xi, p.eta);
	  ref.gradient += coeff * E.Gradient(k, p.xi, p.eta);
	  ref.hessian += coeff * E.Hessian(k, p.xi, p.eta);
	}

      const chemfem::linalg::Matrix2D InvJacT
	= Space->GetMesh().Jacobian(p.cell).Transpose().Invert();

      CachedValues = MapFromReference(ref, InvJacT);
      CachedPoint = p;
      CacheValid = true;

      return CachedValues;
    }

    FEFunction FESpace::Interpolate(ScalarFunction u)
    {
      FEFunction Function(*this);
      Vector Vec(NrDof());

      for(size_t c=0; c<mesh.NrCells(); ++c)
	{
	  const Node& x0 = mesh.Nodes[mesh.Cells[c].LocNode[0]];
	  const chemfem::linalg::Coordinate b{x0.getX(), x0.getY()};

	  const chemfem::linalg::Matrix2D Jac = mesh.Jacobian(c);

	  for(size_t k=0; k<NrLocalDof(); ++k)
	    Vec[GetGlobalIndex(c, k)] = u(b + Jac*refElement.NodalPoint(k));
	}

      Function.SetCoefficients(Vec);

      return Function;
    }

    void FEFunction::WriteVtk(const std::string& filename) const
    {
      Space->GetMesh().WriteVtk(filename, Data);
    }
  }
}
