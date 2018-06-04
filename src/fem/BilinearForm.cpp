#include "fem/BilinearForm.h"
#include "linalg/SparseMatrixInserter.h"
#include "linalg/DenseMatrix.h"
#include "mesh/Mesh.h"
#include "quadrature/QuadFormula.h"

using chemfem::linalg::DenseMatrix;
using chemfem::linalg::SparseMatrix;
using chemfem::linalg::SparseMatrixInserter;

using chemfem::mesh::Mesh;
using chemfem::mesh::Node;
using chemfem::mesh::Cell;

using chemfem::quadrature::QuadratureFormula;
using chemfem::quadrature::QUAD_FORMULA;

namespace chemfem{
  namespace fem{

    double Identity(double x, double y) {return 1.;}
    
    BilinearForm::BilinearForm(const FESpace& TrialSpace, const FESpace& TestSpace)
      : TrialSpace(TrialSpace), TestSpace(TestSpace), Matrix(0,0) {}

    void BilinearForm::AddDiffusionTerm(double (*DiffusionCoeff)(double, double))
    {
      FEExpression expression(SECOND_ORDER, DiffusionCoeff);
      Terms.push_back(expression);
    }

    void BilinearForm::AddLaplaceTerm()
    {

      FEExpression expression(SECOND_ORDER, Identity);
    }

    void BilinearForm::AddReactionTerm(double (*ReactionCoeff)(double, double))
    {
      FEExpression expression(ZERO_ORDER, ReactionCoeff);
      Terms.push_back(expression);
    }

    SparseMatrix& BilinearForm::GetMatrix()
    {
      return Matrix;
    }
        
    void BilinearForm::Assemble()
    {
      Matrix = SparseMatrix(TestSpace.NrDof(), TrialSpace.NrDof());
      SparseMatrixInserter Ins(Matrix);     

      // TODO: Select correct quadrature formula once it is implemented
      QuadratureFormula QuadFormula(QUAD_FORMULA::VERTEX);

      size_t NrQuadPoints = QuadFormula.NrQuadPoints();
      Vector Xi, Eta, Weights;
      QuadFormula.FormulaData(Weights, Xi, Eta);
      
      // Iterate over all cells
      for(std::vector<Cell>::const_iterator cell = TestSpace.mesh.Cells.begin();
	  cell != TestSpace.mesh.Cells.end(); ++cell)
	{
	  double det = cell->Determinant();

	  Node* x0 = cell->LocNode[0];
	  Vector b(2); b[0] = x0->getX(); x0->getY();
	  
	  DenseMatrix Jac = cell->Jacobian();
	  DenseMatrix InvJac(Jac.Invert());

	  // Iterate over all quadrature points
	  Vector::const_iterator Wq, Xiq, Etaq;
	  
	  for(Wq = Weights.begin(), Xiq = Xi.begin(), Etaq = Eta.begin();
	      Wq != Weights.end(); ++Wq, ++Xiq, ++Etaq)
	    {
	      // Determine Quadrature points in world element
	      Vector XiEtaq(2);
	      XiEtaq[0] = *Xiq;
	      XiEtaq[1] = *Etaq;

	      Vector XYq = b + Jac*XiEtaq;
	      	      
	      // Iterate over all terms
	      for(std::vector<FEExpression>::const_iterator Term = Terms.begin();
		  Term != Terms.end(); ++Term)
		{
		  double CoeffVal = Term->Coeff(XYq[0], XYq[1]);

		  
		}
	    }
	}

      Ins.Build();
    }
  }
}
