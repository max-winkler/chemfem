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
      Terms.push_back(expression);
    }

    void BilinearForm::AddReactionTerm(double (*ReactionCoeff)(double, double))
    {
      FEExpression expression(ZERO_ORDER, ReactionCoeff);
      Terms.push_back(expression);
    }

    SparseMatrix& BilinearForm::SystemMatrix()
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

      Vector *GradTest = new Vector[TestSpace.DofPerCell];
      Vector *GradTrial = new Vector[TrialSpace.DofPerCell];
      
      // Iterate over all cells
      int CellInd;
      std::vector<Cell>::const_iterator cell;
      for(cell = TestSpace.mesh.Cells.begin(), CellInd=0;
	  cell != TestSpace.mesh.Cells.end(); ++cell, ++CellInd)
	{
	  double det = cell->Determinant();

	  Node* x0 = cell->LocNode[0];
	  Vector b(2); b[0] = x0->getX(); x0->getY();
	  
	  DenseMatrix Jac = cell->Jacobian();
	  DenseMatrix InvJac(Jac.Invert());

	  DenseMatrix LocMatrix(TestSpace.DofPerCell, TrialSpace.DofPerCell);

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
		  if(&TestSpace.mesh == &TrialSpace.mesh)
		    {
		      double CoeffVal = Term->Coeff(XYq[0], XYq[1]);
		      
		      for(int k=0; k<TestSpace.DofPerCell; ++k)
			GradTest[k] = TestSpace.RefElement().Gradient(k, XYq[0], XYq[1]);
		      for(int l=0; l<TestSpace.DofPerCell; ++l)
			GradTrial[l] = TrialSpace.RefElement().Gradient(l, XYq[0], XYq[1]);
		      
		      switch(Term->Type)
			{
			case SECOND_ORDER:
			  
			  for(int k=0; k<TestSpace.DofPerCell; ++k)
			    for(int l=0; l<TrialSpace.DofPerCell; ++l)
			      LocMatrix[k][l] += 
				(*Wq) * CoeffVal * dot(InvJac*GradTest[k], InvJac*GradTrial[l]) * det;
			  break;
			  /*
			case FIRST_ORDER:

			  break;
			case ZERO_ORDER:

			  break;
			  */
			default:
			  std::cerr << "Assembly of FE expressions of type " << Term->Type
				    << " not implemented yet.\n";
			  return;
			}
		      
		    }
		  else
		    {
		      std::cerr << "Assembly routine for different meshes in trial and test space "
				<< "not implemented yet\n";
		      return;
		    }		  
		} // loop over Terms
	    } // loop over quadrature points

	  // Insert local Matrix into global one
	  for(int k=0; k<TestSpace.DofPerCell; ++k)
	    for(int l=0; l<TestSpace.DofPerCell; ++l)
	      {
		Ins.Insert(TestSpace.GetGlobalIndex(CellInd, k),
			   TrialSpace.GetGlobalIndex(CellInd, l),
			   LocMatrix[k][l]);
	      }
	  
	} // loop over cells

      Ins.Build();

      delete[] GradTest;
      delete[] GradTrial;
    }
  }
}
