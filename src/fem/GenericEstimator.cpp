#include <cmath>
#include <iostream>

#include "fem/GenericEstimator.h"
#include "fem/FESpace.h"

#include "mesh/Mesh.h"
#include "quadrature/QuadFormula.h"

using chemfem::linalg::Coordinate;
using chemfem::linalg::Matrix2D;
using chemfem::linalg::Vector;
using chemfem::linalg::Vector2D;

using chemfem::mesh::Cell;
using chemfem::mesh::CellInfo;
using chemfem::mesh::Edge;
using chemfem::mesh::EdgeType;
using chemfem::mesh::Mesh;
using chemfem::mesh::Node;

using chemfem::quadrature::QuadratureFormula;
using chemfem::quadrature::QUAD_FORMULA;

namespace chemfem{
  namespace fem{

    double Jump(const SolutionState& u, const SolutionState& u_out)
    {
      return u.value - u_out.value;
    }

    double NormalJump(const SolutionState& u, const SolutionState& u_out,
                      const EdgeGeometry& edge)
    {
      return dot(u.gradient - u_out.gradient, edge.normal);
    }

    double GenericEstimator::VolumeResidual::operator()(const Coordinate& pos,
                                                        const CellGeometry& cell,
                                                        const SolutionState& u) const
    {
      // The Laplacian vanishes identically for P1 and is a genuine contribution
      // from P2 on
      const double residual = f(pos) + u.laplacian;

      return cell.h * cell.h * residual * residual;
    }

    double GenericEstimator::EdgeJump::operator()(const Coordinate&,
                                                  const CellGeometry&,
                                                  const EdgeGeometry& edge,
                                                  const SolutionState& u,
                                                  const SolutionState& u_out) const
    {
      const double jump = NormalJump(u, u_out, edge);

      return 0.5 * edge.h * jump * jump;
    }

    GenericEstimator GenericEstimator::Residual(ScalarFunction f)
    {
      GenericEstimator Estimator;

      Estimator.AddVolumeTerm(VolumeResidual(f));
      Estimator.AddEdgeTerm(EdgeJump(), INTERIOR_EDGES);

      return Estimator;
    }

    void GenericEstimator::AddVolumeTerm(VolumeIntegrand term)
    {
      VolumeTerms.push_back(term);
    }

    void GenericEstimator::AddEdgeTerm(EdgeIntegrand term, EdgeSelection selection)
    {
      EdgeTerms.push_back(term);
      EdgeTermSelection.push_back(selection);
    }

    namespace {

      /// State of the discrete solution at a point of the reference element
      void EvalLocal(const FESpace& Space, const FEFunction& u, size_t cell,
                     double xi, double eta, const Matrix2D& InvJacT,
                     SolutionState& state)
      {
        state.value = 0.;
        Vector2D ref_grad;
        Matrix2D ref_hess;

        for(size_t k=0; k<Space.NrLocalDof(); ++k)
          {
            double coeff = u[Space.GetGlobalIndex(cell, k)];

            state.value += coeff * Space.RefElement().Value(k, xi, eta);
            ref_grad += coeff * Space.RefElement().Gradient(k, xi, eta);
            ref_hess += coeff * Space.RefElement().Hessian(k, xi, eta);
          }

        state.gradient = InvJacT * ref_grad;

        // The cells are straight sided, so the reference map x = B*xi + x0 is affine
        // and its Jacobian B is constant. The second derivatives of the map therefore
        // drop out of the chain rule and what remains is the congruence
        //   D^2_x u = B^-T D^2_xi u B^-1,
        // with the very same InvJacT that already transforms the gradient.
        state.hessian = InvJacT * ref_hess * InvJacT.Transpose();
        state.laplacian = state.hessian.Trace();
      }
    }

    Vector GenericEstimator::Assemble(const FEFunction& u) const
    {
      const FESpace& Space = u.GetFESpace();
      const Mesh& mesh = Space.GetMesh();

      const size_t nr_cells = mesh.NrCells();
      Vector Indicators(nr_cells);

      if(VolumeTerms.empty() && EdgeTerms.empty())
        return Indicators;

      QuadratureFormula TriangleQuad(QUAD_FORMULA::GAUSS_7);
      Vector Weights, Xi, Eta;
      TriangleQuad.FormulaData(Weights, Xi, Eta);

      QuadratureFormula LineQuad(QUAD_FORMULA::LINE_GAUSS_5);
      Vector LineWeights, LineNodes, Unused;
      LineQuad.FormulaData(LineWeights, LineNodes, Unused);

      const std::vector<Cell>& Cells = mesh.GetCellList();
      const std::vector<Edge>& Edges = mesh.GetEdgeList();

      for(size_t c=0; c<nr_cells; ++c)
        {
          const Cell& cell = Cells[c];
          CellInfo Info = mesh.GetCellInfo(c);

          const double det = mesh.Determinant(c);
          const double h_T = Info.Diam();

          const Matrix2D Jac = mesh.Jacobian(c);
          const Matrix2D InvJacT = Jac.Transpose().Invert();

          const Node& x0 = mesh.Nodes[cell.LocNode[0]];
          const Coordinate v0{x0.getX(), x0.getY()};

          CellGeometry Geometry;
          Geometry.h = h_T;
          Geometry.area = std::fabs(Info.Volume());
          Geometry.index = c;

          double value = 0.;

          // ---------------------------------------------------------------- volume
          for(size_t q=0; q<Weights.size(); ++q)
            {
              const Coordinate RefPoint{Xi[q], Eta[q]};
              const Coordinate pos = v0 + Jac*RefPoint;

              SolutionState State;
              EvalLocal(Space, u, c, Xi[q], Eta[q], InvJacT, State);

              for(size_t t=0; t<VolumeTerms.size(); ++t)
                value += Weights[q] * VolumeTerms[t](pos, Geometry, State)
                       * std::fabs(det);
            }

          // ------------------------------------------------------------------ edges
          if(!EdgeTerms.empty())
            {
              for(int k=0; k<3; ++k)
                {
                  const Edge& edge = Edges[cell.LocEdge[k]];
                  const bool boundary = (edge.Type() != EdgeType::INTERFACE_EDGE);

                  size_t neigh_cell = c;
                  int neigh_edge = k;

                  if(!boundary)
                    {
                      neigh_cell = edge.GetNeighbor(c);
                      neigh_edge = Cells[neigh_cell].EdgeIndex(cell.LocEdge[k]);
                    }

                  const Matrix2D NeighJac = mesh.Jacobian(neigh_cell);
                  const Matrix2D NeighInvJacT = NeighJac.Transpose().Invert();

                  EdgeGeometry EdgeGeom;
                  EdgeGeom.h = Info.EdgeLength(k);
                  EdgeGeom.normal = Info.Normal(k);
                  EdgeGeom.local_index = k;
                  EdgeGeom.boundary = boundary;

                  for(size_t q=0; q<LineWeights.size(); ++q)
                    {
                      const double s = LineNodes[q];

                      double xi, eta;
                      EdgeToRefCoords(k, s, xi, eta);

                      const Coordinate RefPoint{xi, eta};
                      const Coordinate pos = v0 + Jac*RefPoint;

                      SolutionState State, NeighState;
                      EvalLocal(Space, u, c, xi, eta, InvJacT, State);

                      if(boundary)
                        NeighState = State;
                      else
                        {
                          // The neighbor traverses the shared edge in the opposite
                          // direction, so the same physical point sits at 1-s there.
                          double neigh_xi, neigh_eta;
                          EdgeToRefCoords(neigh_edge, 1.-s, neigh_xi, neigh_eta);

                          EvalLocal(Space, u, neigh_cell, neigh_xi, neigh_eta,
                                    NeighInvJacT, NeighState);
                        }

                      for(size_t t=0; t<EdgeTerms.size(); ++t)
                        {
                          const EdgeSelection selection = EdgeTermSelection[t];

                          if(boundary && selection == INTERIOR_EDGES) continue;
                          if(!boundary && selection == BOUNDARY_EDGES) continue;

                          value += LineWeights[q]
                                 * EdgeTerms[t](pos, Geometry, EdgeGeom, State, NeighState)
                                 * EdgeGeom.h;
                        }
                    }
                }
            }

          Indicators[c] = value;
        }

      return Indicators;
    }
  }
}
