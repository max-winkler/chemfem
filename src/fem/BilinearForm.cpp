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
using chemfem::mesh::CellInfo;
using chemfem::mesh::EdgeType;

using chemfem::quadrature::QuadratureFormula;
using chemfem::quadrature::QUAD_FORMULA;

namespace chemfem{
  namespace fem{

    double Identity(const Coordinate&) {return 1.;}

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

    void BilinearForm::AddVolumeTerm(Integrand term)
    {
      VolumeTerms.push_back(term);
    }

    void BilinearForm::AddVolumeTerm(PointIntegrand term)
    {
      PointVolumeTerms.push_back(term);
    }

    void BilinearForm::AddBoundaryTerm(Integrand term, BoundaryIndicator part)
    {
      BoundaryTerms.push_back(BoundaryTerm{term, nullptr, part});
    }

    void BilinearForm::AddBoundaryTerm(BoundaryIntegrand term, BoundaryIndicator part)
    {
      BoundaryTerms.push_back(BoundaryTerm{nullptr, term, part});
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

    void BilinearForm::Assemble(SparseMatrixInserter& Ins, size_t RowOffset, size_t ColOffset,
                                bool Transposed)
    {
      Assemble(Ins, std::vector<Placement>(1, Placement{RowOffset, ColOffset, Transposed}));
    }

    void BilinearForm::Assemble(SparseMatrixInserter& Ins,
                                const std::vector<Placement>& Places)
    {
      if(&TestSpace.mesh != &TrialSpace.mesh)
        {
          std::cerr << "Assembly routine for different meshes in trial and test space "
                    << "not implemented yet\n";
          return;
        }

      const Mesh& mesh = TestSpace.mesh;

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

      const bool HasIntegrands = !VolumeTerms.empty() || !PointVolumeTerms.empty();

      // The basis functions on the reference element in the quadrature points, the same
      // for every cell
      std::vector<PointValues> RefTest, RefTrial;
      if(HasIntegrands)
        {
          RefTest = TabulateReference(TestSpace.RefElement(), Xi, Eta);
          RefTrial = TabulateReference(TrialSpace.RefElement(), Xi, Eta);
        }

      std::vector<PointValues> TestValues(NrTest), TrialValues(NrTrial);

      // Iterate over all cells
      int CellInd;
      std::vector<Cell>::const_iterator cell;
      for(cell = mesh.Cells.begin(), CellInd=0;
          cell != mesh.Cells.end(); ++cell, ++CellInd)
        {
          double det = mesh.Determinant(CellInd);

          const Node& x0 = mesh.Nodes[cell->LocNode[0]];
          const Coordinate b{x0.getX(), x0.getY()};

          const Matrix2D Jac = mesh.Jacobian(CellInd);
          const Matrix2D InvJac = Jac.Transpose().Invert();

          DenseMatrix LocMatrix(NrTest, NrTrial);

          // Iterate over all quadrature points
          Vector::const_iterator Wq, Xiq, Etaq;
          size_t q;

          for(Wq = Weights.begin(), Xiq = Xi.begin(), Etaq = Eta.begin(), q = 0;
              Wq != Weights.end(); ++Wq, ++Xiq, ++Etaq, ++q)
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

                    default:
                      std::cerr << "Assembly of FE expressions of type " << Term->Type
                                << " not implemented yet.\n";
                      return;
                    }
                } // loop over Terms

              if(HasIntegrands)
                {
                  for(int k=0; k<NrTest; ++k)
                    TestValues[k] = MapFromReference(RefTest[q*NrTest + k], InvJac);

                  for(int l=0; l<NrTrial; ++l)
                    TrialValues[l] = MapFromReference(RefTrial[q*NrTrial + l], InvJac);

                  for(size_t t=0; t<VolumeTerms.size(); ++t)
                    for(int k=0; k<NrTest; ++k)
                      for(int l=0; l<NrTrial; ++l)
                        LocMatrix[k][l] += (*Wq)
                          * VolumeTerms[t](TrialValues[l], TestValues[k]) * det;

                  if(!PointVolumeTerms.empty())
                    {
                      const QuadPoint Point{XYq, size_t(CellInd), *Xiq, *Etaq};

                      for(size_t t=0; t<PointVolumeTerms.size(); ++t)
                        for(int k=0; k<NrTest; ++k)
                          for(int l=0; l<NrTrial; ++l)
                            LocMatrix[k][l] += (*Wq)
                              * PointVolumeTerms[t](Point, TrialValues[l], TestValues[k]) * det;
                    }
                }
            } // loop over quadrature points

          InsertLocalMatrix(Ins, Places, CellInd, LocMatrix);

        } // loop over cells

      if(!BoundaryTerms.empty())
        {
          QuadratureFormula LineFormula(QUAD_FORMULA::LINE_GAUSS_5);

          Vector LineWeights, LineNodes, Unused;
          LineFormula.FormulaData(LineWeights, LineNodes, Unused);

          for(size_t e=0; e<mesh.Edges.size(); ++e)
            {
              if(mesh.Edges[e].Type() != EdgeType::BOUNDARY_EDGE)
                continue;

              const size_t CellIndex = mesh.Edges[e].GetNeighbor(-1);
              const Cell& EdgeCell = mesh.Cells[CellIndex];
              const int LocEdge = EdgeCell.EdgeIndex(e);

              const Node& P0 = mesh.Nodes[EdgeCell.LocNode[LocEdge]];
              const Node& P1 = mesh.Nodes[EdgeCell.LocNode[(LocEdge+1)%3]];
              const Coordinate Midpoint{0.5*(P0.getX() + P1.getX()),
                                        0.5*(P0.getY() + P1.getY())};

              std::vector<size_t> Active;
              for(size_t t=0; t<BoundaryTerms.size(); ++t)
                if(!BoundaryTerms[t].part || BoundaryTerms[t].part(Midpoint))
                  Active.push_back(t);

              if(Active.empty())
                continue;

              CellInfo Info = mesh.GetCellInfo(CellIndex);

              EdgeGeometry EdgeGeom;
              EdgeGeom.h = Info.EdgeLength(LocEdge);
              EdgeGeom.normal = Info.Normal(LocEdge);
              EdgeGeom.local_index = LocEdge;
              EdgeGeom.boundary = true;

              const Matrix2D InvJac = mesh.Jacobian(CellIndex).Transpose().Invert();

              DenseMatrix LocMatrix(NrTest, NrTrial);

              for(size_t q=0; q<LineWeights.size(); ++q)
                {
                  const double s = LineNodes[q];

                  double xi, eta;
                  EdgeToRefCoords(LocEdge, s, xi, eta);

                  const QuadPoint Point{Coordinate{(1.-s)*P0.getX() + s*P1.getX(),
                                                   (1.-s)*P0.getY() + s*P1.getY()},
                                        CellIndex, xi, eta};

                  for(int k=0; k<NrTest; ++k)
                    TestValues[k] = MapFromReference(
                      ReferenceValues(TestSpace.RefElement(), k, xi, eta), InvJac);

                  for(int l=0; l<NrTrial; ++l)
                    TrialValues[l] = MapFromReference(
                      ReferenceValues(TrialSpace.RefElement(), l, xi, eta), InvJac);

                  for(size_t a=0; a<Active.size(); ++a)
                    {
                      const BoundaryTerm& T = BoundaryTerms[Active[a]];

                      for(int k=0; k<NrTest; ++k)
                        for(int l=0; l<NrTrial; ++l)
                          {
                            const double value = T.integrand
                              ? T.integrand(TrialValues[l], TestValues[k])
                              : T.point_integrand(Point, EdgeGeom, TrialValues[l], TestValues[k]);

                            LocMatrix[k][l] += LineWeights[q] * value * EdgeGeom.h;
                          }
                    }
                }

              InsertLocalMatrix(Ins, Places, CellIndex, LocMatrix);
            }
        }

      delete[] GradTest;
      delete[] GradTrial;
      delete[] ValueTest;
      delete[] ValueTrial;
    }

    void BilinearForm::InsertLocalMatrix(SparseMatrixInserter& Ins,
                                         const std::vector<Placement>& Places, size_t CellInd,
                                         const DenseMatrix& LocMatrix)
    {
      const DofManager& TestDofs = TestSpace.Dofs;
      const DofManager& TrialDofs = TrialSpace.Dofs;

      const int NrTest = TestSpace.NrLocalDof();
      const int NrTrial = TrialSpace.NrLocalDof();

      for(int k=0; k<NrTest; ++k)
        for(int l=0; l<NrTrial; ++l)
          {
            size_t DofTrial = TrialDofs.GlobalIndex(CellInd, l);
            size_t DofTest = TestDofs.GlobalIndex(CellInd, k);

            // Skip for test functions not in the test
            if(!TestDofs.IsFree(DofTest))
              continue;

            if(TrialDofs.IsFree(DofTrial))
              {
                // DOF is a free DOF
                const size_t Row = TestDofs.ReducedIndex(DofTest);
                const size_t Col = TrialDofs.ReducedIndex(DofTrial);

                for(size_t p=0; p<Places.size(); ++p)
                  {
                    if(Places[p].transposed)
                      Ins.Insert(Places[p].row + Col, Places[p].col + Row, LocMatrix[k][l]);
                    else
                      Ins.Insert(Places[p].row + Row, Places[p].col + Col, LocMatrix[k][l]);
                  }
              }
            else
              {
                //DOF is a Dirichlet DOF
                // \todo Modify this when implementing inhomogeneous Dirichlet conditions
                double Value = 0.;
                DirichletRhs[TestDofs.ReducedIndex(DofTest)] += LocMatrix[k][l] * Value;
              }
          }
    }
  }
}
