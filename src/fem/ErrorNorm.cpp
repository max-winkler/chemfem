#include <cmath>

#include "fem/ErrorNorm.h"
#include "quadrature/QuadFormula.h"

using chemfem::linalg::DenseMatrix;
using chemfem::mesh::Node;
using chemfem::mesh::Cell;
using chemfem::quadrature::QuadratureFormula;

namespace chemfem{
  namespace fem{

    void ErrorNorm::SetExactValue(ScalarFunction* Value)
    {
      this->Value = Value;
    }

    void ErrorNorm::SetExactGradient(VectorFunction* Gradient)
    {
      this->Gradient = Gradient;
    }

    void ErrorNorm::SetFEFunction(FEFunction& FESolution)
    {
      this->FESolution = &FESolution;
    }

    double ErrorNorm::Error(Norm norm) const
    {
      double error = 0.;

      const FESpace& Space = FESolution->GetFESpace();
      const Mesh& mesh = Space.GetMesh();

      QuadratureFormula quad(chemfem::quadrature::GAUSS_7);
      Vector Weights, Xi, Eta;
      quad.FormulaData(Weights, Xi, Eta);
      
      std::vector<Cell>::const_iterator it_cell;
      for(it_cell = mesh.GetCellList().begin();
	  it_cell != mesh.GetCellList().end(); ++it_cell)
	{
	  const Cell& cell = *it_cell;

	  double det = cell.Determinant();

	  const Node& x0 = cell.GetNode(0);
	  Vector b(2); b[0] = x0.getX(); x0.getY();
	  
	  DenseMatrix Jac = cell.Jacobian();

	  double loc_error = 0.;

	  const size_t* LocalDof = Space.GetLocalDofMap(cell.Index());
	  
	  Vector::const_iterator Xiq, Etaq, Wq;
	  for(Xiq = Xi.begin(), Etaq = Eta.begin(), Wq = Weights.begin();
	      Xiq != Xi.end(); ++Xiq, ++Etaq, ++Wq)
	    {
	      Vector XiEtaq(2);
	      XiEtaq[0] = *Xiq, XiEtaq[1] = *Etaq;

	      Vector XYq = b + Jac*XiEtaq;
	      
	      if(norm == L2 || norm == H1)
		{
		  // Value of FE solution
		  double fe_value = 0.;
		  for(int k=0; k<Space.RefElement().NrDof(); ++k)
		    fe_value += (*FESolution)[LocalDof[k]] * Space.RefElement().Value(k, *Xiq, *Etaq);

		  // Value of exact solution
		  double ex_value = (*Value)(XYq[0], XYq[1]);
		  loc_error += 0.5 * (*Wq) * pow(fe_value - ex_value, 2.);		  
		}
	    } // loop over quadrature points
	  
	  error += loc_error * det;
	    
	} // loop over cells
      
      error = sqrt(error);
      
      return error;
    }
    
  }
}
