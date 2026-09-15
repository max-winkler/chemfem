#include <cmath>

#include "fem/ErrorNorm.h"
#include "quadrature/QuadFormula.h"

using chemfem::linalg::DenseMatrix;
using chemfem::mesh::Mesh;
using chemfem::mesh::Node;
using chemfem::mesh::Cell;
using chemfem::quadrature::QuadratureFormula;

namespace chemfem{
  namespace fem{

    ErrorNorm::ErrorNorm(ScalarFunction Value) : Value(Value) {}

    ErrorNorm::ErrorNorm(ScalarFunction Value, VectorFunction Gradient)
      : Value(Value), Gradient(Gradient) {}

    void ErrorNorm::SetExactValue(ScalarFunction Value)
    {
      this->Value = Value;
    }

    void ErrorNorm::SetExactGradient(VectorFunction Gradient)
    {
      this->Gradient = Gradient;
    }

    void ErrorNorm::SetFEFunction(FEFunction& FESolution)
    {
      this->FESolution = &FESolution;
    }

    double ErrorNorm::Compute(Norm norm) const
    {
      return Compute(*FESolution, norm);
    }

    double ErrorNorm::Compute(const FEFunction& Solution, Norm norm, int component) const
    {
      double error = 0.;

      const FESpace& Space = Solution.GetFESpace();
      const Mesh& mesh = Space.GetMesh();

      // A Piola mapped basis function contributes to both components, and it has no
      // gradient, so only the L2 norm is available for such a space
      const bool Piola = Space.RefElement().Mapping() == ContravariantPiola;

      if(Piola && norm != L2)
        {
          std::cerr << "Error: Only the L2 norm is available on a space with a Piola "
                    << "mapping.\n";
          return 0.;
        }

      QuadratureFormula quad(chemfem::quadrature::GAUSS_7);
      Vector Weights, Xi, Eta;
      quad.FormulaData(Weights, Xi, Eta);
      
      std::vector<Cell>::const_iterator it_cell;
      size_t CellInd;
        for(it_cell = mesh.GetCellList().begin(), CellInd = 0;
	  it_cell != mesh.GetCellList().end(); ++it_cell, ++CellInd)
	{
	  const Cell& cell = *it_cell;

	  const Node& x0 = mesh.Nodes[cell.LocNode[0]];
	  const chemfem::linalg::Coordinate b{x0.getX(), x0.getY()};

	  const chemfem::linalg::Matrix2D Jac = mesh.Jacobian(CellInd);
	  const chemfem::linalg::Matrix2D InvJac = Jac.Transpose().Invert();
	  
	  double det = Jac.Determinant();
	  
	  double loc_error = 0.;
	  
	  const size_t* LocalDof = Space.GetLocalDofMap(CellInd);
	  	  
	  Vector::const_iterator Xiq, Etaq, Wq;
	  for(Xiq = Xi.begin(), Etaq = Eta.begin(), Wq = Weights.begin();
	      Xiq != Xi.end(); ++Xiq, ++Etaq, ++Wq)
	    {
	      const chemfem::linalg::Coordinate XiEtaq{*Xiq, *Etaq};
	      const chemfem::linalg::Coordinate XYq = b + Jac*XiEtaq;
	      
	      if(norm == L2 || norm == H1)
		{
		  // Value of FE solution
		  double fe_value = 0.;
		  for(int k=0; k<Space.RefElement().NrDof(); ++k)
		    {
		      // On a vector valued space only the basis functions of the component
		      // that is asked for contribute
		      if(!Piola && Space.RefElement().Component(k) != component)
			continue;

		      double form_value = Piola
			? MapFromReference(Space.RefElement().VectorReference(k, *Xiq, *Etaq),
					   Jac, det, Space.LocalSign(CellInd, k)).value[component]
			: Space.RefElement().Value(k, *Xiq, *Etaq);
		      double dof_value = Solution[LocalDof[k]];
		      fe_value += dof_value * form_value;
		    }
		  // Value of exact solution
		  double ex_value = Value(XYq);

		  double diff = fe_value - ex_value;
		  loc_error += (*Wq) * pow(diff, 2.);
		}
	      if(norm == H1 || norm == H1_SEMI)
		{
		  chemfem::linalg::Vector2D fe_grad;
		  for(int k=0; k<Space.RefElement().NrDof(); ++k)
		    {
		      // On a vector valued space only the basis functions of the component
		      // that is asked for contribute
		      if(!Piola && Space.RefElement().Component(k) != component)
			continue;

		      double dof_value = Solution[LocalDof[k]];
		      fe_grad += dof_value * Space.RefElement().Gradient(k, *Xiq, *Etaq);
		    }

		  chemfem::linalg::Vector2D ex_grad = Gradient(XYq);

		  chemfem::linalg::Vector2D diff = ex_grad - InvJac*fe_grad;
		  loc_error += (*Wq) * dot(diff, diff);
		}
	    } // loop over quadrature points
	  
	  error += loc_error * det;
	    
	} // loop over cells
      
      error = sqrt(error);
      
      return error;
    }
    
  }
}
