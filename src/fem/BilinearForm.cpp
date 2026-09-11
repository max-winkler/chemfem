#include "fem/BilinearForm.h"
#include "linalg/SparseMatrixInserter.h"
#include "linalg/DenseMatrix.h"
#include "mesh/Mesh.h"
#include "quadrature/QuadFormula.h"

using chemfem::linalg::Coordinate;
using chemfem::linalg::DenseMatrix;
using chemfem::linalg::Matrix2D;
using chemfem::linalg::Vector2D;
using chemfem::linalg::SparseMatrix;
using chemfem::linalg::SparseMatrixInserter;

using chemfem::mesh::Mesh;
using chemfem::mesh::Node;
using chemfem::mesh::Cell;

using chemfem::quadrature::QuadratureFormula;
using chemfem::quadrature::QUAD_FORMULA;

namespace chemfem{
  namespace fem{

    double Identity(const Coordinate&) {return 1.;}

    namespace {

      double Apply(FEOperator op, double value, const Vector2D& grad)
      {
        switch(op)
          {
          case DX: return grad.x;
          case DY: return grad.y;
          default: return value;
          }
      }
    }

    BilinearForm::BilinearForm(const FESpace& TrialSpace, const FESpace& TestSpace)
      : TrialSpace(TrialSpace), TestSpace(TestSpace), Matrix(0,0),
	DirichletRhs(TestSpace.NrDof()) {}

    void BilinearForm::AddDiffusionTerm(ScalarFunction DiffusionCoeff)
    {
      FEExpression expression(SECOND_ORDER, DiffusionCoeff);
      Terms.push_back(expression);
    }

    void BilinearForm::AddLaplaceTerm()
    {
      FEExpression expression(SECOND_ORDER, Identity);
      Terms.push_back(expression);
    }

    void BilinearForm::AddConvectionTerm(VectorFunction ConvectionField)
    {
      FEExpression expression(FIRST_ORDER, ConvectionField);
      Terms.push_back(expression);
    }

    void BilinearForm::AddReactionTerm(ScalarFunction ReactionCoeff)
    {
      FEExpression expression(ZERO_ORDER, ReactionCoeff);
      Terms.push_back(expression);
    }

    void BilinearForm::AddTerm(ScalarFunction Coeff, FEOperator TrialOp, FEOperator TestOp)
    {
      FEExpression expression(Coeff, TrialOp, TestOp);
      Terms.push_back(expression);
    }

    void BilinearForm::AddTerm(FEOperator TrialOp, FEOperator TestOp)
    {
      AddTerm(Identity, TrialOp, TestOp);
    }

    SparseMatrix& BilinearForm::SystemMatrix()
    {
      return Matrix;
    }

    const FESpace& BilinearForm::GetTrialSpace() const
    {
      return TrialSpace;
    }

    const FESpace& BilinearForm::GetTestSpace() const
    {
      return TestSpace;
    }

    void BilinearForm::Assemble()
    {
      Matrix = SparseMatrix(TestSpace.NrFreeDof(), TrialSpace.NrFreeDof());
      SparseMatrixInserter Ins(Matrix);

      Assemble(Ins, 0, 0);

      Ins.Build();
    }

    void BilinearForm::Assemble(SparseMatrixInserter& Ins, size_t RowOffset, size_t ColOffset)
    {
      const DofManager& TestDofs = TestSpace.Dofs;
      const DofManager& TrialDofs = TrialSpace.Dofs;

      const int NrTest = TestSpace.NrLocalDof();
      const int NrTrial = TrialSpace.NrLocalDof();

      // TODO: Select correct quadrature formula once it is implemented
      QuadratureFormula QuadFormula(QUAD_FORMULA::GAUSS_7);

      Vector Xi, Eta, Weights;
      QuadFormula.FormulaData(Weights, Xi, Eta);

      Vector2D *GradTest = new Vector2D[NrTest];
      Vector2D *GradTrial = new Vector2D[NrTrial];

      double *ValueTest = new double[NrTest];
      double *ValueTrial = new double[NrTrial];

      // Iterate over all cells
      int CellInd;
      std::vector<Cell>::const_iterator cell;
      for(cell = TestSpace.mesh.Cells.begin(), CellInd=0;
          cell != TestSpace.mesh.Cells.end(); ++cell, ++CellInd)
        {
          double det = TestSpace.mesh.Determinant(CellInd);

          Node& x0 = TestSpace.mesh.Nodes[cell->LocNode[0]];
          const Coordinate b{x0.getX(), x0.getY()};

          const Matrix2D Jac = TestSpace.mesh.Jacobian(CellInd);
          const Matrix2D InvJac = Jac.Transpose().Invert();

          DenseMatrix LocMatrix(NrTest, NrTrial);

          // Iterate over all quadrature points
          Vector::const_iterator Wq, Xiq, Etaq;

          for(Wq = Weights.begin(), Xiq = Xi.begin(), Etaq = Eta.begin();
              Wq != Weights.end(); ++Wq, ++Xiq, ++Etaq)
            {
              // Determine Quadrature points in world element
              const Coordinate XiEtaq{*Xiq, *Etaq};
              const Coordinate XYq = b + Jac*XiEtaq;

              for(int k=0; k<NrTest; ++k)
                {
                  GradTest[k] = InvJac * TestSpace.RefElement().Gradient(k, *Xiq, *Etaq);
                  ValueTest[k] = TestSpace.RefElement().Value(k, *Xiq, *Etaq);
                }

              for(int l=0; l<NrTrial; ++l)
                {
                  GradTrial[l] = InvJac * TrialSpace.RefElement().Gradient(l, *Xiq, *Etaq);
                  ValueTrial[l] = TrialSpace.RefElement().Value(l, *Xiq, *Etaq);
                }

              // Iterate over all terms
              for(std::vector<FEExpression>::const_iterator Term = Terms.begin();
                  Term != Terms.end(); ++Term)
                {
                  if(&TestSpace.mesh == &TrialSpace.mesh)
                    {
                      switch(Term->Type)
                        {
                        case SECOND_ORDER:
                          {
                            const double CoeffVal = Term->Coeff(XYq);

                            for(int k=0; k<NrTest; ++k)
                              for(int l=0; l<NrTrial; ++l)
                                {
                                  LocMatrix[k][l] +=
                                    (*Wq) * CoeffVal
                                    * dot(GradTest[k], GradTrial[l])
                                    * det;
                                }
                          }
                          break;

                        case FIRST_ORDER:
                          {
                            const Vector2D ConvectionField = Term->VecCoeff(XYq);

                            for(int k=0; k<NrTest; ++k)
                              for(int l=0; l<NrTrial; ++l)
                                {
                                  LocMatrix[k][l] +=
                                    (*Wq) * dot(ConvectionField, GradTrial[l]) * ValueTest[k]
                                    * det;
                                }
                          }
                          break;

                        case ZERO_ORDER:
                          {
                            const double CoeffVal = Term->Coeff(XYq);

                            for(int k=0; k<NrTest; ++k)
                              for(int l=0; l<NrTrial; ++l)
                                LocMatrix[k][l] += (*Wq) * CoeffVal * ValueTest[k]
                                  * ValueTrial[l] * det;
                          }
                          break;

                        case GENERAL:
                          {
                            const double CoeffVal = Term->Coeff(XYq);

                            for(int k=0; k<NrTest; ++k)
                              for(int l=0; l<NrTrial; ++l)
                                LocMatrix[k][l] += (*Wq) * CoeffVal
                                  * Apply(Term->TrialOp, ValueTrial[l], GradTrial[l])
                                  * Apply(Term->TestOp, ValueTest[k], GradTest[k]) * det;
                          }
                          break;

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
          for(int k=0; k<NrTest; ++k)
            for(int l=0; l<NrTrial; ++l)
              {
                size_t DofTrial = TrialDofs.GlobalIndex(CellInd, l);
                size_t DofTest = TestDofs.GlobalIndex(CellInd, k);

                // Skip for test functions not in the test
                if(!TestDofs.IsFree(DofTest))
                  continue;

                if(TrialDofs.IsFree(DofTrial))
                  // DOF is a free DOF
                  Ins.Insert(RowOffset + TestDofs.ReducedIndex(DofTest),
                             ColOffset + TrialDofs.ReducedIndex(DofTrial), LocMatrix[k][l]);
                else
                  {
                    //DOF is a Dirichlet DOF
                    // \todo Modify this when implementing inhomogeneous Dirichlet conditions
                    double Value = 0.;
                    DirichletRhs[TestDofs.ReducedIndex(DofTest)] += LocMatrix[k][l] * Value;
                  }
              }

        } // loop over cells

      delete[] GradTest;
      delete[] GradTrial;
      delete[] ValueTest;
      delete[] ValueTrial;
    }
  }
}
