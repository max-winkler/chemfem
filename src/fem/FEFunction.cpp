#include <cmath>
#include <vector>

#include "fem/FEFunction.h"
#include "fem/DirichletValues.h"
#include "fem/DofTransform.h"

#include "quadrature/QuadFormula.h"

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

    const Vector& FEFunction::Coefficients() const
    {
      return Data;
    }

    void FEFunction::CreateFunction(const Vector& FreeDof)
    {
      Data = Space->IncorporateBC(FreeDof);
      CacheValid = NOTHING;
    }

    void FEFunction::CreateFunction(const Vector& FreeDof, const DirichletValues& g)
    {
      Data = Space->IncorporateBC(FreeDof, g);
      CacheValid = NOTHING;
    }

    void FEFunction::SetCoefficients(const Vector& Data)
    {
      this->Data = Data;
      CacheValid = NOTHING;
    }

    namespace {

      bool SamePoint(const QuadPoint& a, const QuadPoint& b)
      {
	return a.cell == b.cell && a.xi == b.xi && a.eta == b.eta;
      }
    }

    double FEFunction::Value(const QuadPoint& p) const
    {
      if((CacheValid == VALUE_ONLY || CacheValid == EVERYTHING) && SamePoint(p, CachedPoint))
	return CachedValues.value;

      const ScalarElement* E = Space->AsScalar();

      if(!E)
	{
	  std::cerr << "Error: This space has vector valued shape functions, use "
		    << "VectorValue.\n";
	  return 0.;
	}

      std::vector<double> Coeff;
      GatherLocalCoefficients(*Space, Data, p.cell, Coeff);

      double value = 0.;
      for(size_t k=0; k<Space->NrLocalDof(); ++k)
	value += Coeff[k] * E->Value(k, p.xi, p.eta);

      CachedValues.value = value;
      CachedPoint = p;
      CacheValid = VALUE_ONLY;

      return value;
    }

    chemfem::linalg::Vector2D FEFunction::VectorValue(const QuadPoint& p) const
    {
      if(CacheValid == VECTOR_ONLY && SamePoint(p, CachedPoint))
	return CachedVector;

      const Element& E = Space->RefElement();

      chemfem::linalg::Vector2D value;

      std::vector<double> Coeff;
      GatherLocalCoefficients(*Space, Data, p.cell, Coeff);

      if(const VectorElement* Vec = Space->AsVector())
	{
	  const chemfem::linalg::Matrix2D Jac = Space->GetMesh().Jacobian(p.cell);
	  const double det = Space->GetMesh().Determinant(p.cell);

	  for(size_t k=0; k<Space->NrLocalDof(); ++k)
	    {
	      const chemfem::linalg::Vector2D basis
		= MapFromReference(ReferenceVector{Vec->Value(k, p.xi, p.eta),
						   Vec->Gradient(k, p.xi, p.eta)},
				   Jac, det, Space->LocalSign(p.cell, k)).value;

	      const double coeff = Coeff[k];

	      value[0] += coeff * basis[0];
	      value[1] += coeff * basis[1];
	    }
	}
      else
	for(size_t k=0; k<Space->NrLocalDof(); ++k)
	  value[E.Component(k)] += Coeff[k]
	    * Space->AsScalar()->Value(k, p.xi, p.eta);

      CachedVector = value;
      CachedPoint = p;
      CacheValid = VECTOR_ONLY;

      return value;
    }

    PointValues FEFunction::Evaluate(const QuadPoint& p) const
    {
      if(CacheValid == EVERYTHING && SamePoint(p, CachedPoint))
	return CachedValues;

      const ScalarElement* E = Space->AsScalar();

      if(!E)
	{
	  std::cerr << "Error: Evaluate needs scalar shape functions.\n";
	  return PointValues{0., chemfem::linalg::Vector2D(), chemfem::linalg::Matrix2D(), 0.};
	}

      std::vector<double> Coeff;
      GatherLocalCoefficients(*Space, Data, p.cell, Coeff);

      PointValues ref;
      ref.value = 0.;

      for(size_t k=0; k<Space->NrLocalDof(); ++k)
	{
	  const double coeff = Coeff[k];

	  ref.value += coeff * E->Value(k, p.xi, p.eta);
	  ref.gradient += coeff * E->Gradient(k, p.xi, p.eta);
	  ref.hessian += coeff * E->Hessian(k, p.xi, p.eta);
	}

      const chemfem::linalg::Matrix2D InvJacT
	= Space->GetMesh().Jacobian(p.cell).Transpose().Invert();

      CachedValues = MapFromReference(ref, InvJacT);
      CachedPoint = p;
      CacheValid = EVERYTHING;

      return CachedValues;
    }

    double FEFunction::Mean() const
    {
      const chemfem::mesh::Mesh& mesh = Space->GetMesh();

      chemfem::quadrature::QuadratureFormula Quad(chemfem::quadrature::QUAD_FORMULA::GAUSS_7);

      Vector Weights, Xi, Eta;
      Quad.FormulaData(Weights, Xi, Eta);

      double integral = 0., area = 0.;

      std::vector<double> Coeff;

      for(size_t c=0; c<mesh.NrCells(); ++c)
	{
	  const double det = std::fabs(mesh.Determinant(c));

	  GatherLocalCoefficients(*Space, Data, c, Coeff);

	  for(size_t q=0; q<Weights.size(); ++q)
	    {
	      double value = 0.;
	      for(size_t k=0; k<Space->NrLocalDof(); ++k)
		value += Coeff[k]
		  * Space->AsScalar()->Value(k, Xi[q], Eta[q]);

	      integral += Weights[q] * value * det;
	      area += Weights[q] * det;
	    }
	}

      return integral/area;
    }

    void FEFunction::SubtractMean()
    {
      const double mean = Mean();

      for(size_t k=0; k<Data.size(); ++k)
	Data[k] -= mean;

      CacheValid = NOTHING;
    }

    FEFunction FESpace::Interpolate(ScalarFunction u)
    {
      FEFunction Function(*this);
      Vector Vec(NrDof());

      for(int k=0; k<refElement.NrDof(); ++k)
	if(refElement.Dof(k).type != DofType::PointValue)
	  {
	    std::cerr << "Error: Interpolation sets every DOF to a function value, which is "
		      << "wrong for an element whose DOFs are derivatives or edge moments. "
		      << "The zero function is returned.\n";
	    return Function;
	  }

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
      if(Space->RefElement().DofsPerVertex() != 1)
	{
	  std::cerr << "Error: This writes the first coefficient of every node, which is the "
		    << "value in that vertex only for an element with one DOF per vertex. Use "
		    << "VtkOutput for this space.\n";
	  return;
	}

      Space->GetMesh().WriteVtk(filename, Data);
    }
  }
}
