#include "fem/LinearForm.h"

#include "quadrature/QuadFormula.h"

using chemfem::linalg::Vector;
using chemfem::linalg::Vector2D;
using chemfem::linalg::Matrix2D;
using chemfem::linalg::DenseMatrix;
using chemfem::linalg::Coordinate;

using chemfem::mesh::Mesh;
using chemfem::mesh::Node;
using chemfem::mesh::Cell;
using chemfem::mesh::CellInfo;
using chemfem::mesh::Edge;
using chemfem::mesh::EdgeType;

using chemfem::quadrature::QuadratureFormula;
using chemfem::quadrature::QUAD_FORMULA;

namespace chemfem{
  namespace fem{

    LinearForm::LinearForm(const FESpace& TestSpace) : TestSpace(TestSpace) {}

    void LinearForm::AddVolumeForce(ScalarFunction F)
    {
      FEExpression Expression(VOLUME_FORCE, F);
      Terms.push_back(Expression);
    }

    void LinearForm::AddNeumannBC(ScalarFunction G)
    {
      FEExpression Expression(NEUMANN_BC, G);
      Terms.push_back(Expression);
    }

    void LinearForm::AddVolumeTerm(Integrand term)
    {
      VolumeTerms.push_back(term);
    }

    void LinearForm::AddVolumeTerm(PointIntegrand term)
    {
      PointVolumeTerms.push_back(term);
    }

    void LinearForm::AddBoundaryTerm(Integrand term, BoundaryIndicator part)
    {
      BoundaryTerms.push_back(BoundaryTerm{term, nullptr, part});
    }

    void LinearForm::AddBoundaryTerm(BoundaryIntegrand term, BoundaryIndicator part)
    {
      BoundaryTerms.push_back(BoundaryTerm{nullptr, term, part});
    }

    Vector& LinearForm::LoadVector()
    {
      return Vec;
    }

    const FESpace& LinearForm::GetTestSpace() const
    {
      return TestSpace;
    }

    void LinearForm::Assemble()
    {
      Vec = Vector(TestSpace.NrFreeDof());

      const Mesh& mesh = TestSpace.mesh;
      const DofManager& Dofs = TestSpace.Dofs;
      const int NrTest = TestSpace.NrLocalDof();

      // TODO: Select correct quadrature formula once it is implemented
      QuadratureFormula QuadFormula(QUAD_FORMULA::GAUSS_7);

      Vector Xi, Eta, Weights;
      QuadFormula.FormulaData(Weights, Xi, Eta);

      double *TestFuncValue = new double[NrTest];

      const bool HasIntegrands = !VolumeTerms.empty() || !PointVolumeTerms.empty();

      // The basis functions on the reference element in the quadrature points, the same
      // for every cell
      std::vector<PointValues> RefTest;
      if(HasIntegrands)
	RefTest = TabulateReference(TestSpace.RefElement(), Xi, Eta);

      std::vector<PointValues> TestValues(NrTest);

      // Iterate over all cells
      int CellInd;
      std::vector<Cell>::const_iterator cell;
      for(cell = mesh.Cells.begin(), CellInd=0;
	  cell != mesh.Cells.end(); ++cell, ++CellInd)
	{
	  double det = mesh.Determinant(CellInd);

	  const Node& x0 = mesh.Nodes[cell->LocNode[0]];
	  const chemfem::linalg::Coordinate b{x0.getX(), x0.getY()};

	  const chemfem::linalg::Matrix2D Jac = mesh.Jacobian(CellInd);
	  const Matrix2D InvJac = Jac.Transpose().Invert();

	  Vector LocVec(NrTest);

	  // Iterate over all quadrature points
	  Vector::const_iterator Wq, Xiq, Etaq;
	  size_t q;

	  for(Wq = Weights.begin(), Xiq = Xi.begin(), Etaq = Eta.begin(), q = 0;
	      Wq != Weights.end(); ++Wq, ++Xiq, ++Etaq, ++q)
	    {
	      // Determine Quadrature points in world element
	      const chemfem::linalg::Coordinate XiEtaq{*Xiq, *Etaq};
	      const chemfem::linalg::Coordinate XYq = b + Jac*XiEtaq;

	      // Function value of test functions
	      for(int k=0; k<NrTest; ++k)
		TestFuncValue[k] = TestSpace.RefElement().Value(k, *Xiq, *Etaq);

	      // Iterate over all terms
	      for(std::vector<FEExpression>::const_iterator Term = Terms.begin();
		  Term != Terms.end(); ++Term)
		{
		  switch(Term->Type)
		    {
		    case VOLUME_FORCE:
		      {
			const double CoeffVal = Term->Coeff(XYq);

			for(int k=0; k<NrTest; ++k)
			  LocVec[k] += (*Wq) * CoeffVal * TestFuncValue[k] * det;
		      }
		      break;

		    case NEUMANN_BC:
		      // Assembled below by the loop over the boundary edges
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

		  for(size_t t=0; t<VolumeTerms.size(); ++t)
		    for(int k=0; k<NrTest; ++k)
		      LocVec[k] += (*Wq) * VolumeTerms[t](TestValues[k]) * det;

		  if(!PointVolumeTerms.empty())
		    {
		      const QuadPoint Point{XYq, size_t(CellInd), *Xiq, *Etaq};

		      for(size_t t=0; t<PointVolumeTerms.size(); ++t)
			for(int k=0; k<NrTest; ++k)
			  LocVec[k] += (*Wq) * PointVolumeTerms[t](Point, TestValues[k]) * det;
		    }
		}
	    } // loop over quadrature points
	  for(int k=0; k<NrTest; ++k)
	    {
	      size_t GlobalIndex = TestSpace.GetGlobalIndex(CellInd, k);
	      if(Dofs.IsFree(GlobalIndex))
		Vec[Dofs.ReducedIndex(GlobalIndex)] += LocVec[k];
	    }

	} // loop over cells

      delete[] TestFuncValue;

      // Neumann boundary conditions on the boundary edges that are not Dirichlet edges
      QuadratureFormula LineFormula(QUAD_FORMULA::LINE_GAUSS_5);

      Vector LineWeights, LineNodes, Unused;
      LineFormula.FormulaData(LineWeights, LineNodes, Unused);

      const std::vector<Edge>& Edges = mesh.Edges;

      for(size_t e=0; e<Edges.size(); ++e)
	{
	  if(Edges[e].Type() != EdgeType::BOUNDARY_EDGE || Dofs.IsDirichletEdge(e))
	    continue;

	  const size_t CellIndex = Edges[e].GetNeighbor(-1);
	  const Cell& EdgeCell = mesh.Cells[CellIndex];
	  const int k = EdgeCell.EdgeIndex(e);

	  const Node& P0 = mesh.Nodes[EdgeCell.LocNode[k]];
	  const Node& P1 = mesh.Nodes[EdgeCell.LocNode[(k+1)%3]];
	  const double length = P0.Dist(P1);

	  Vector LocVec(NrTest);

	  for(size_t q=0; q<LineWeights.size(); ++q)
	    {
	      const double s = LineNodes[q];
	      const Coordinate XYq{(1.-s)*P0.getX() + s*P1.getX(),
				   (1.-s)*P0.getY() + s*P1.getY()};

	      double xi, eta;
	      EdgeToRefCoords(k, s, xi, eta);

	      for(std::vector<FEExpression>::const_iterator Term = Terms.begin();
		  Term != Terms.end(); ++Term)
		{
		  if(Term->Type != NEUMANN_BC)
		    continue;

		  const double CoeffVal = Term->Coeff(XYq);

		  for(int i=0; i<NrTest; ++i)
		    LocVec[i] += LineWeights[q] * CoeffVal
		      * TestSpace.RefElement().Value(i, xi, eta) * length;
		}
	    }

	  for(int i=0; i<NrTest; ++i)
	    {
	      size_t GlobalIndex = TestSpace.GetGlobalIndex(CellIndex, i);
	      if(Dofs.IsFree(GlobalIndex))
		Vec[Dofs.ReducedIndex(GlobalIndex)] += LocVec[i];
	    }
	}

      // Boundary terms of AddBoundaryTerm, on the parts of the boundary they are given for
      for(size_t e=0; e<Edges.size() && !BoundaryTerms.empty(); ++e)
	{
	  if(Edges[e].Type() != EdgeType::BOUNDARY_EDGE)
	    continue;

	  const size_t CellIndex = Edges[e].GetNeighbor(-1);
	  const Cell& EdgeCell = mesh.Cells[CellIndex];
	  const int LocEdge = EdgeCell.EdgeIndex(e);

	  const Node& P0 = mesh.Nodes[EdgeCell.LocNode[LocEdge]];
	  const Node& P1 = mesh.Nodes[EdgeCell.LocNode[(LocEdge+1)%3]];
	  const Coordinate Midpoint{0.5*(P0.getX() + P1.getX()), 0.5*(P0.getY() + P1.getY())};

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

	  Vector LocVec(NrTest);

	  for(size_t q=0; q<LineWeights.size(); ++q)
	    {
	      const double s = LineNodes[q];

	      double xi, eta;
	      EdgeToRefCoords(LocEdge, s, xi, eta);

	      const QuadPoint Point{Coordinate{(1.-s)*P0.getX() + s*P1.getX(),
					       (1.-s)*P0.getY() + s*P1.getY()},
				    CellIndex, xi, eta};

	      for(int i=0; i<NrTest; ++i)
		TestValues[i] = MapFromReference(
		  ReferenceValues(TestSpace.RefElement(), i, xi, eta), InvJac);

	      for(size_t a=0; a<Active.size(); ++a)
		{
		  const BoundaryTerm& T = BoundaryTerms[Active[a]];

		  for(int i=0; i<NrTest; ++i)
		    {
		      const double value = T.integrand
			? T.integrand(TestValues[i])
			: T.point_integrand(Point, EdgeGeom, TestValues[i]);

		      LocVec[i] += LineWeights[q] * value * EdgeGeom.h;
		    }
		}
	    }

	  for(int i=0; i<NrTest; ++i)
	    {
	      size_t GlobalIndex = TestSpace.GetGlobalIndex(CellIndex, i);
	      if(Dofs.IsFree(GlobalIndex))
		Vec[Dofs.ReducedIndex(GlobalIndex)] += LocVec[i];
	    }
	}
    }


  }
}
