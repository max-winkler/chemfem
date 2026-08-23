#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

#include "fem/GenericEstimator.h"

#include "fem/LinearForm.h"
#include "fem/BilinearForm.h"
#include "fem/LagrangeElement.h"
#include "fem/FEFunction.h"
#include "mesh/UnitSquareMesh.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// ---------------------------------------------------------------------------------
// Model problem
//
//   -Laplace(u) = f   in  Omega = (0,1)^2,      u = 0  on  dOmega
//
// with the exact solution u = 16 x(1-x) y(1-y). For this operator the energy norm
// is the H1 seminorm, which is what the residual estimator below bounds.
// ---------------------------------------------------------------------------------
double f(const Coordinate& p)
{
  return 32.*(p.x*(1.-p.x) + p.y*(1.-p.y));
}

Vector2D exact_grad(const Coordinate& p)
{
  Vector2D grad;
  grad[0] = 16.*p.y*(1.-p.y)*(1.-2.*p.x);
  grad[1] = 16.*p.x*(1.-p.x)*(1.-2.*p.y);
  return grad;
}

// ---------------------------------------------------------------------------------
// The standard residual based a posteriori error estimator for the H1 seminorm
//
//   eta_T^2 = h_T^2 ||f + Laplace(u_h)||_{L2(T)}^2
//           + 1/2 sum_{E in dT, E interior} h_E ||[du_h/dn]_E||_{L2(E)}^2
//
//   eta^2 = sum_T eta_T^2
//
// The factor 1/2 distributes each interior edge over the two cells sharing it: an
// interior edge is visited once from either side.
//
// Reliability   |u - u_h|_H1 <= C_rel  eta
// Efficiency    eta          <= C_eff (|u - u_h|_H1 + oscillation)
//
// so eta/|u - u_h|_H1 has to settle at a constant, which is what the table below
// checks. The terms are handed out separately as well, to see how the two
// contributions balance.
// ---------------------------------------------------------------------------------

/// h_T^2 |f + Laplace(u_h)|^2
struct VolumeResidual
{
  double operator()(const Coordinate& pos, const CellGeometry& cell,
                    const SolutionState& u) const
  {
    // For P1 the Laplacian of u_h vanishes on every cell, u.laplacian is 0 then.
    double residual = f(pos) + u.laplacian;

    return cell.h * cell.h * residual * residual;
  }
};

/// 1/2 h_E |[du_h/dn]|^2
struct EdgeJump
{
  double operator()(const Coordinate& pos, const CellGeometry& cell,
                    const EdgeGeometry& edge,
                    const SolutionState& u, const SolutionState& u_out) const
  {
    double jump = NormalJump(u, u_out, edge);

    return 0.5 * edge.h * jump * jump;
  }
};

/// |u - u_h|_{H1(T)}^2, the exact error the estimator is measured against
struct ExactErrorH1
{
  double operator()(const Coordinate& pos, const CellGeometry& cell,
                    const SolutionState& u) const
  {
    Vector2D grad = exact_grad(pos);

    double dx = u.gradient[0] - grad[0];
    double dy = u.gradient[1] - grad[1];

    return dx*dx + dy*dy;
  }
};

// --------------------------------------------------------------------- self tests
struct One
{
  double operator()(const Coordinate&, const CellGeometry&,
                    const SolutionState&) const { return 1.; }
  double operator()(const Coordinate&, const CellGeometry&,
                    const EdgeGeometry&,
                    const SolutionState&, const SolutionState&) const { return 1.; }
};

struct MeshSize
{
  double operator()(const Coordinate& pos, const CellGeometry& cell,
                    const SolutionState& u) const { return cell.h; }
};

struct ValueJump
{
  double operator()(const Coordinate& pos, const CellGeometry& cell,
                    const EdgeGeometry& edge,
                    const SolutionState& u, const SolutionState& u_out) const { return std::fabs(Jump(u, u_out)); }
};

// ---------------------------------------------------------------------------------

double Sum(const Vector& v)
{
  double s = 0.;
  for(Vector::const_iterator it = v.begin(); it != v.end(); ++it)
    s += *it;
  return s;
}

FEFunction Solve(FESpace& Space)
{
  BilinearForm A(Space, Space);
  A.AddLaplaceTerm();
  A.Assemble();

  LinearForm F(Space);
  F.AddVolumeForce(f);
  F.Assemble();

  Vector X(A.SystemMatrix().Solve(F.LoadVector()));

  FEFunction Sol(Space);
  Sol.CreateFunction(X);

  return Sol;
}

/// One row of the tables. Collected during the loops, because the CG solver writes
/// its own progress to stdout.
struct Row
{
  size_t cells;
  double eta, eta_volume, eta_edge, error;
};

void PrintTable(const std::string& caption, const std::vector<Row>& rows)
{
  std::cout << "\n" << caption << "\n";
  std::cout << std::string(78, '=') << "\n";
  std::cout << std::setw(8)  << "Zellen"
            << std::setw(12) << "eta"
            << std::setw(12) << "eta_Vol"
            << std::setw(12) << "eta_Kante"
            << std::setw(12) << "|u-u_h|_H1"
            << std::setw(12) << "eta/Fehler"
            << std::setw(10) << "eoc" << std::endl;

  std::cout << std::setprecision(4) << std::fixed;

  for(size_t i=0; i<rows.size(); ++i)
    {
      std::cout << std::setw(8)  << rows[i].cells
                << std::setw(12) << rows[i].eta
                << std::setw(12) << rows[i].eta_volume
                << std::setw(12) << rows[i].eta_edge
                << std::setw(12) << rows[i].error
                << std::setw(12) << rows[i].eta / rows[i].error;

      if(i > 0)
        {
          // Rate with respect to the number of cells. Halving h quadruples that
          // number, so this is the usual order of convergence in h.
          double ratio = double(rows[i].cells) / double(rows[i-1].cells);
          std::cout << std::setw(10) << 2.*log(rows[i-1].eta/rows[i].eta)/log(ratio);
        }

      std::cout << std::endl;
    }

  std::cout << std::defaultfloat << std::setprecision(6);
}

int main()
{
  // ------------------------------------------------------------------------------
  // 0) Does the quadrature integrate what it claims to integrate?
  //
  //    A term returning 1 has to reproduce the area of the domain resp. the length
  //    of its boundary. On the criss-cross mesh of the unit square every triangle
  //    has diameter 1/(n-1), so the integral of h_T is 1/(n-1) as well. The last
  //    check verifies that both sides of an interior edge are parametrized
  //    consistently: u_h is continuous, so its jump has to vanish.
  // ------------------------------------------------------------------------------
  {
    const size_t n = 5;

    UnitSquareMesh mesh(n);
    LagrangeElement element(1);
    FESpace Space(mesh, element);
    FEFunction Zero(Space);

    GenericEstimator Area;
    Area.AddVolumeTerm(One());

    GenericEstimator Perimeter;
    Perimeter.AddEdgeTerm(One(), BOUNDARY_EDGES);

    GenericEstimator Size;
    Size.AddVolumeTerm(MeshSize());

    GenericEstimator Continuity;
    Continuity.AddEdgeTerm(ValueJump(), INTERIOR_EDGES);

    FEFunction Sol = Solve(Space);

    std::cout << "0) Quadratur\n" << std::string(78, '=') << "\n";
    std::cout << "   Flaeche des Einheitsquadrats        : " << Sum(Area.Assemble(Zero))
              << "   (exakt 1)\n";
    std::cout << "   Laenge des Randes                   : " << Sum(Perimeter.Assemble(Zero))
              << "   (exakt 4)\n";
    std::cout << "   Integral von h_T                    : " << Sum(Size.Assemble(Zero))
              << "   (exakt " << 1./(n-1) << ")\n";
    std::cout << "   Sprung von u_h ueber inneren Kanten : " << Sum(Continuity.Assemble(Sol))
              << "   (exakt 0)\n";
  }

  // ------------------------------------------------------------------------------
  // The estimator is built once. It stores only the terms, not a solution, so the
  // same object is evaluated for every mesh of both loops below.
  // ------------------------------------------------------------------------------
  GenericEstimator Estimator;
  Estimator.AddVolumeTerm(VolumeResidual());
  Estimator.AddEdgeTerm(EdgeJump(), INTERIOR_EDGES);

  GenericEstimator VolumePart;
  VolumePart.AddVolumeTerm(VolumeResidual());

  GenericEstimator EdgePart;
  EdgePart.AddEdgeTerm(EdgeJump(), INTERIOR_EDGES);

  GenericEstimator ExactError;
  ExactError.AddVolumeTerm(ExactErrorH1());

  std::vector<Row> uniform, adaptive;

  // ------------------------------------------------------------------------------
  // 1) Uniform refinement. Both eta and the exact error have to converge with order
  //    one, and their ratio has to settle at a constant.
  // ------------------------------------------------------------------------------
  {
    UnitSquareMesh mesh(2);

    for(int level=0; level<6; ++level)
      {
        LagrangeElement element(1);
        FESpace Space(mesh, element);

        FEFunction Sol = Solve(Space);

        Row row;
        row.cells      = mesh.NrCells();
        row.eta        = sqrt(Sum(Estimator.Assemble(Sol)));
        row.eta_volume = sqrt(Sum(VolumePart.Assemble(Sol)));
        row.eta_edge   = sqrt(Sum(EdgePart.Assemble(Sol)));
        row.error      = sqrt(Sum(ExactError.Assemble(Sol)));
        uniform.push_back(row);

        mesh.RefineUniform();
        mesh.RefineUniform();
      }
  }

  // ------------------------------------------------------------------------------
  // 2) The same estimator driving an adaptive loop.
  // ------------------------------------------------------------------------------
  {
    UnitSquareMesh mesh(3);

    for(int level=0; level<8; ++level)
      {
        LagrangeElement element(1);
        FESpace Space(mesh, element);

        FEFunction Sol = Solve(Space);

        Vector Indicators = Estimator.Assemble(Sol);

        Row row;
        row.cells      = mesh.NrCells();
        row.eta        = sqrt(Sum(Indicators));
        row.eta_volume = sqrt(Sum(VolumePart.Assemble(Sol)));
        row.eta_edge   = sqrt(Sum(EdgePart.Assemble(Sol)));
        row.error      = sqrt(Sum(ExactError.Assemble(Sol)));
        adaptive.push_back(row);

        // Maximum strategy, as in tests/AdaptivityTest.cpp
        double MaxIndicator = 0.;
        for(Vector::const_iterator it = Indicators.begin(); it != Indicators.end(); ++it)
          if(*it > MaxIndicator) MaxIndicator = *it;

        std::vector<bool> Marker(mesh.NrCells(), false);
        for(size_t c=0; c<mesh.NrCells(); ++c)
          Marker[c] = (Indicators[c] > 0.5*MaxIndicator);

        mesh.Refine(Marker);

        if(!mesh.Check())
          {
            std::cerr << "ERROR: mesh is broken after refinement step " << level+1 << ".\n";
            return 1;
          }
      }
  }

  PrintTable("1) Uniforme Verfeinerung", uniform);
  PrintTable("2) Adaptive Verfeinerung, gesteuert von denselben Indikatoren", adaptive);

  std::cout << "\nEstimatorTest was successful.\n";
  return 0;
}
