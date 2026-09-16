#include "fem/BilinearForm.h"

#include "fem/DofTransform.h"
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
	DirichletTerm(TestSpace.NrFreeDof()) {}

    void BilinearForm::SetDirichletValues(const DirichletValues& g)
    {
      if(&g.GetFESpace() != &TrialSpace)
        {
          std::cerr << "Error: The Dirichlet values have to belong to the trial space of the "
                    << "bilinear form.\n";
          return;
        }

      PrescribedValues = &g;
    }

    const Vector& BilinearForm::DirichletRhs() const
    {
      return DirichletTerm;
    }

    namespace {

      bool CheckVector(const FESpace& Space, const char* which)
      {
        if(Space.IsVectorValued())
          return true;

        std::cerr << "Error: This term needs a vector valued " << which
                  << " space, but that space is scalar.\n";
        return false;
      }

      bool CheckScalar(const FESpace& Space, const char* which)
      {
        if(!Space.IsVectorValued())
          return true;

        std::cerr << "Error: This term needs a scalar " << which
                  << " space, but that space is vector valued.\n";
        return false;
      }

      /// Both spaces of a term that is evaluated with scalar values
      bool CheckScalarPair(const FESpace& Trial, const FESpace& Test)
      {
        return CheckScalar(Trial, "trial") && CheckScalar(Test, "test");
      }
    }

    void BilinearForm::AddDiffusionTerm(ScalarFunction DiffusionCoeff)
    {
      if(!CheckScalarPair(TrialSpace, TestSpace))
        return;

      FEExpression expression(SECOND_ORDER, DiffusionCoeff);
      Terms.push_back(expression);
    }

    void BilinearForm::AddLaplaceTerm()
    {
      if(!CheckScalarPair(TrialSpace, TestSpace))
        return;

      FEExpression expression(SECOND_ORDER, Identity);
      Terms.push_back(expression);
    }

    void BilinearForm::AddConvectionTerm(VectorFunction ConvectionField)
    {
      if(!CheckScalarPair(TrialSpace, TestSpace))
        return;

      FEExpression expression(FIRST_ORDER, ConvectionField);
      Terms.push_back(expression);
    }

    void BilinearForm::AddReactionTerm(ScalarFunction ReactionCoeff)
    {
      if(!CheckScalarPair(TrialSpace, TestSpace))
        return;

      FEExpression expression(ZERO_ORDER, ReactionCoeff);
      Terms.push_back(expression);
    }

    void BilinearForm::AddVolumeTerm(Integrand term)
    {
      if(CheckScalarPair(TrialSpace, TestSpace))
        VolumeTerms.push_back(term);
    }

    void BilinearForm::AddVolumeTerm(PointIntegrand term)
    {
      if(CheckScalarPair(TrialSpace, TestSpace))
        PointVolumeTerms.push_back(term);
    }

    void BilinearForm::AddVolumeTerm(VectorIntegrand term)
    {
      if(CheckVector(TrialSpace, "trial") && CheckVector(TestSpace, "test"))
        VectorTerms.push_back(term);
    }

    void BilinearForm::AddVolumeTerm(VectorPointIntegrand term)
    {
      if(CheckVector(TrialSpace, "trial") && CheckVector(TestSpace, "test"))
        VectorPointTerms.push_back(term);
    }

    void BilinearForm::AddVolumeTerm(VectorScalarIntegrand term)
    {
      if(CheckVector(TrialSpace, "trial") && CheckScalar(TestSpace, "test"))
        VectorScalarTerms.push_back(term);
    }

    void BilinearForm::AddVolumeTerm(ScalarVectorIntegrand term)
    {
      if(CheckScalar(TrialSpace, "trial") && CheckVector(TestSpace, "test"))
        ScalarVectorTerms.push_back(term);
    }

    void BilinearForm::AddBoundaryTerm(Integrand term, BoundaryIndicator part)
    {
      // The boundary loop builds scalar values only, it would ignore the components
      if(CheckScalarPair(TrialSpace, TestSpace))
        BoundaryTerms.push_back(BoundaryTerm{term, nullptr, part});
    }

    void BilinearForm::AddBoundaryTerm(BoundaryIntegrand term, BoundaryIndicator part)
    {
      if(CheckScalarPair(TrialSpace, TestSpace))
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
      DirichletTerm = Vector(TestSpace.NrFreeDof());

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

      const bool ScalarTest = !VolumeTerms.empty() || !PointVolumeTerms.empty()
        || !VectorScalarTerms.empty();
      const bool VectorTest = !VectorTerms.empty() || !VectorPointTerms.empty()
        || !ScalarVectorTerms.empty();
      const bool ScalarTrial = !VolumeTerms.empty() || !PointVolumeTerms.empty()
        || !ScalarVectorTerms.empty();
      const bool VectorTrial = !VectorTerms.empty() || !VectorPointTerms.empty()
        || !VectorScalarTerms.empty();

      const bool HasIntegrands = ScalarTest || VectorTest;
      const bool NeedsPoint = !PointVolumeTerms.empty() || !VectorPointTerms.empty();

      const VectorElement* VecTest = TestSpace.AsVector();
      const VectorElement* VecTrial = TrialSpace.AsVector();

      const bool PiolaTest = VecTest != nullptr;
      const bool PiolaTrial = VecTrial != nullptr;

      // The basis functions on the reference element in the quadrature points, the same
      // for every cell
      std::vector<PointValues> RefTest, RefTrial;
      std::vector<ReferenceVector> VecRefTest, VecRefTrial;
      if(HasIntegrands)
        {
          if(PiolaTest)
            VecRefTest = TabulateVectorReference(*VecTest, Xi, Eta);
          else
            RefTest = TabulateReference(*TestSpace.AsScalar(), Xi, Eta);

          if(PiolaTrial)
            VecRefTrial = TabulateVectorReference(*VecTrial, Xi, Eta);
          else
            RefTrial = TabulateReference(*TrialSpace.AsScalar(), Xi, Eta);
        }

      std::vector<PointValues> TestValues(NrTest), TrialValues(NrTrial);
      std::vector<VectorValues> TestVectors(NrTest), TrialVectors(NrTrial);

      // Asked once, like the tabulation above, so the cell loop does not walk the DOFs of
      // every element again only to learn that there is nothing to transform
      const bool TransformTest = NeedsDofTransform(TestSpace);
      const bool TransformTrial = NeedsDofTransform(TrialSpace);

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

              // Only the FEExpression terms read these, and an element with a Piola mapping
              // has no scalar values to offer
              if(!Terms.empty())
                {
                  for(int k=0; k<NrTest; ++k)
                    {
                      GradTest[k] = InvJac * TestSpace.AsScalar()->Gradient(k, *Xiq, *Etaq);
                      ValueTest[k] = TestSpace.AsScalar()->Value(k, *Xiq, *Etaq);
                    }

                  for(int l=0; l<NrTrial; ++l)
                    {
                      GradTrial[l] = InvJac * TrialSpace.AsScalar()->Gradient(l, *Xiq, *Etaq);
                      ValueTrial[l] = TrialSpace.AsScalar()->Value(l, *Xiq, *Etaq);
                    }
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
                    {
                      if(PiolaTest)
                        {
                          TestVectors[k] = MapFromReference(VecRefTest[q*NrTest + k], Jac, det,
                                                            TestSpace.LocalSign(CellInd, k));
                          continue;
                        }

                      const PointValues& ref = RefTest[q*NrTest + k];

                      if(ScalarTest)
                        TestValues[k] = MapFromReference(ref, InvJac);
                      if(VectorTest)
                        TestVectors[k] = MapFromReference(ref, InvJac,
                                                          TestSpace.RefElement().Component(k));
                    }

                  for(int l=0; l<NrTrial; ++l)
                    {
                      if(PiolaTrial)
                        {
                          TrialVectors[l] = MapFromReference(VecRefTrial[q*NrTrial + l], Jac, det,
                                                             TrialSpace.LocalSign(CellInd, l));
                          continue;
                        }

                      const PointValues& ref = RefTrial[q*NrTrial + l];

                      if(ScalarTrial)
                        TrialValues[l] = MapFromReference(ref, InvJac);
                      if(VectorTrial)
                        TrialVectors[l] = MapFromReference(ref, InvJac,
                                                           TrialSpace.RefElement().Component(l));
                    }

                  for(size_t t=0; t<VolumeTerms.size(); ++t)
                    for(int k=0; k<NrTest; ++k)
                      for(int l=0; l<NrTrial; ++l)
                        LocMatrix[k][l] += (*Wq)
                          * VolumeTerms[t](TrialValues[l], TestValues[k]) * det;

                  for(size_t t=0; t<VectorTerms.size(); ++t)
                    for(int k=0; k<NrTest; ++k)
                      for(int l=0; l<NrTrial; ++l)
                        LocMatrix[k][l] += (*Wq)
                          * VectorTerms[t](TrialVectors[l], TestVectors[k]) * det;

                  for(size_t t=0; t<VectorScalarTerms.size(); ++t)
                    for(int k=0; k<NrTest; ++k)
                      for(int l=0; l<NrTrial; ++l)
                        LocMatrix[k][l] += (*Wq)
                          * VectorScalarTerms[t](TrialVectors[l], TestValues[k]) * det;

                  for(size_t t=0; t<ScalarVectorTerms.size(); ++t)
                    for(int k=0; k<NrTest; ++k)
                      for(int l=0; l<NrTrial; ++l)
                        LocMatrix[k][l] += (*Wq)
                          * ScalarVectorTerms[t](TrialValues[l], TestVectors[k]) * det;

                  if(NeedsPoint)
                    {
                      const QuadPoint Point{XYq, size_t(CellInd), *Xiq, *Etaq};

                      for(size_t t=0; t<PointVolumeTerms.size(); ++t)
                        for(int k=0; k<NrTest; ++k)
                          for(int l=0; l<NrTrial; ++l)
                            LocMatrix[k][l] += (*Wq)
                              * PointVolumeTerms[t](Point, TrialValues[l], TestValues[k]) * det;

                      for(size_t t=0; t<VectorPointTerms.size(); ++t)
                        for(int k=0; k<NrTest; ++k)
                          for(int l=0; l<NrTrial; ++l)
                            LocMatrix[k][l] += (*Wq)
                              * VectorPointTerms[t](Point, TrialVectors[l], TestVectors[k]) * det;
                    }
                }
            } // loop over quadrature points

          if(TransformTest)
            TransformLocalRows(LocMatrix, NrTrial, TestSpace, CellInd, Jac);
          if(TransformTrial)
            TransformLocalColumns(LocMatrix, NrTest, TrialSpace, CellInd, Jac);

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
                      ReferenceValues(*TestSpace.AsScalar(), k, xi, eta), InvJac);

                  for(int l=0; l<NrTrial; ++l)
                    TrialValues[l] = MapFromReference(
                      ReferenceValues(*TrialSpace.AsScalar(), l, xi, eta), InvJac);

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

              const Matrix2D CellJac = mesh.Jacobian(CellIndex);

              if(TransformTest)
                TransformLocalRows(LocMatrix, NrTrial, TestSpace, CellIndex, CellJac);
              if(TransformTrial)
                TransformLocalColumns(LocMatrix, NrTest, TrialSpace, CellIndex, CellJac);

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
            else if(PrescribedValues)
              {
                // The DOF is prescribed, so its column moves to the right hand side
                DirichletTerm[TestDofs.ReducedIndex(DofTest)]
                  += LocMatrix[k][l] * (*PrescribedValues)[DofTrial];
              }
          }
    }
  }
}
