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
      double error = 0.;

      const FESpace& Space = FESolution->GetFESpace();
      const Mesh& mesh = Space.GetMesh();

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
	  Vector b(2);
	  b[0] = x0.getX(); b[1] = x0.getY();
	  
	  DenseMatrix Jac = mesh.Jacobian(CellInd);
	  DenseMatrix InvJac(Jac.Transpose().Invert());
	  
	  double det = Jac.Determinant();
	  
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
		    {
		      double form_value = Space.RefElement().Value(k, *Xiq, *Etaq);
		      double dof_value = (*FESolution)[LocalDof[k]];
		      fe_value += dof_value * form_value;
		    }
		  // Value of exact solution
		  double ex_value = Value(XYq[0], XYq[1]);

		  double diff = fe_value - ex_value;
		  loc_error += (*Wq) * pow(diff, 2.);
		}
	      if(norm == H1 || norm == H1_SEMI)
		{
		  Vector fe_grad(2);
		  for(int k=0; k<Space.RefElement().NrDof(); ++k)
		    {
		      Vector form_grad = Space.RefElement().Gradient(k, *Xiq, *Etaq);
		      double dof_value = (*FESolution)[LocalDof[k]];
		      fe_grad += dof_value*form_grad;
		    }
		  
		  Vector ex_grad = Gradient(XYq[0], XYq[1]);
		  
		  Vector diff = ex_grad - InvJac*fe_grad;
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
