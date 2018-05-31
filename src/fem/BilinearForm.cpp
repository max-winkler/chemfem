#include "fem/BilinearForm.h"
#include "linalg/SparseMatrixInserter.h"
#include "linalg/DenseMatrix.h"
#include "mesh/Mesh.h"

using chemfem::linalg::DenseMatrix;
using chemfem::linalg::SparseMatrix;
using chemfem::linalg::SparseMatrixInserter;

using chemfem::mesh::Mesh;
using chemfem::mesh::Cell;

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

      // Iterate over all cells
      for(std::vector<Cell>::const_iterator cell = TestSpace.mesh.Cells.begin();
	  cell != TestSpace.mesh.Cells.end(); ++cell)
	{
	  double det = cell->Determinant();
	  DenseMatrix Jac = cell->Jacobian();
	  DenseMatrix InvJac(Jac.Invert());

	  
	}

      Ins.Build();
    }
  }
}
