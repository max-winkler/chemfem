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
// Model problem with a variable diffusion coefficient
//
//   -div(a grad u) = f   in  Omega = (0,1)^2,      u = 0  on  dOmega
//
//   a(x,y) = 1 + contrast * x * y
//   u(x,y) = 16 x(1-x) y(1-y)
//
// Two things are demonstrated. First, the contrast is a member of the functor and
// therefore a runtime value - which is what the former raw function pointer could
// not express. Second, every function of the coordinates is written on Coordinate,
// so the same object feeds the bilinear form, the right hand side and the error
// estimator without any adapter in between.
// ---------------------------------------------------------------------------------

/// The diffusion coefficient a, carrying its own parameter
struct Diffusion
{
  double contrast;

  explicit Diffusion(double contrast) : contrast(contrast) {}

  double operator()(const Coordinate& p) const
  {
    return 1. + contrast*p.x*p.y;
  }

  Vector2D Gradient(const Coordinate& p) const
  {
    Vector2D g;
    g[0] = contrast*p.y;
    g[1] = contrast*p.x;
    return g;
  }
};

Vector2D grad_u_exact(const Coordinate& p)
{
  Vector2D grad;
  grad[0] = 16.*p.y*(1.-p.y)*(1.-2.*p.x);
  grad[1] = 16.*p.x*(1.-p.x)*(1.-2.*p.y);
  return grad;
}

double laplace_u_exact(const Coordinate& p)
{
  return -32.*(p.x*(1.-p.x) + p.y*(1.-p.y));
}

/// f = -div(a grad u) = -(grad a . grad u) - a Laplace(u)
struct RightHandSide
{
  Diffusion a;

  explicit RightHandSide(const Diffusion& a) : a(a) {}

  double operator()(const Coordinate& p) const
  {
    Vector2D ga = a.Gradient(p);
    Vector2D gu = grad_u_exact(p);

    return -(ga[0]*gu[0] + ga[1]*gu[1]) - a(p)*laplace_u_exact(p);
  }
};

// ---------------------------------------------------------------------------------
// The weighted residual estimator for the diffusion problem
//
//   eta_T^2 = h_T^2/a ||f + div(a grad u_h)||_{L2(T)}^2
//           + 1/2 sum_E h_E/a ||[a du_h/dn]_E||_{L2(E)}^2
//
// measured against the energy norm |||v|||^2 = int a |grad v|^2. For P1 the
// Laplacian of u_h vanishes, so div(a grad u_h) reduces to grad a . grad u_h. The
// coefficient is continuous here, so the weight is the same on both sides of an
// edge.
// ---------------------------------------------------------------------------------
struct WeightedVolumeResidual
{
  Diffusion a;
  RightHandSide f;

  explicit WeightedVolumeResidual(const Diffusion& a) : a(a), f(a) {}

  double operator()(const Coordinate& pos, const CellGeometry& cell,
                    const SolutionState& u) const
  {
    Vector2D ga = a.Gradient(pos);

    double div_flux = ga[0]*u.gradient[0] + ga[1]*u.gradient[1]
                    + a(pos)*u.laplacian;
    double residual = f(pos) + div_flux;

    return cell.h * cell.h / a(pos) * residual * residual;
  }
};

struct WeightedEdgeJump
{
  Diffusion a;

  explicit WeightedEdgeJump(const Diffusion& a) : a(a) {}

  double operator()(const Coordinate& pos, const CellGeometry& cell,
                    const EdgeGeometry& edge,
                    const SolutionState& u, const SolutionState& u_out) const
  {
    double jump = a(pos) * NormalJump(u, u_out, edge);

    return 0.5 * edge.h / a(pos) * jump * jump;
  }
};

/// |||u - u_h|||^2 on the cell
struct EnergyError
{
  Diffusion a;

  explicit EnergyError(const Diffusion& a) : a(a) {}

  double operator()(const Coordinate& pos, const CellGeometry& cell,
                    const SolutionState& u) const
  {
    Vector2D gu = grad_u_exact(pos);

    double dx = u.gradient[0] - gu[0];
    double dy = u.gradient[1] - gu[1];

    return a(pos) * (dx*dx + dy*dy);
  }
};

// ---------------------------------------------------------------------------------

double Sum(const Vector& v)
{
  double s = 0.;
  for(Vector::const_iterator it = v.begin(); it != v.end(); ++it)
    s += *it;
  return s;
}

struct Row
{
  size_t cells;
  double eta, error;
};

void PrintTable(const std::string& caption, const std::vector<Row>& rows)
{
  std::cout << "\n" << caption << "\n" << std::string(60, '=') << "\n";
  std::cout << std::setw(8)  << "Zellen"
            << std::setw(14) << "eta"
            << std::setw(14) << "|||u-u_h|||"
            << std::setw(14) << "eta/Fehler"
            << std::setw(10) << "eoc" << std::endl;

  std::cout << std::setprecision(4) << std::fixed;

  for(size_t i=0; i<rows.size(); ++i)
    {
      std::cout << std::setw(8)  << rows[i].cells
                << std::setw(14) << rows[i].eta
                << std::setw(14) << rows[i].error
                << std::setw(14) << rows[i].eta / rows[i].error;

      if(i > 0)
        {
          double ratio = double(rows[i].cells) / double(rows[i-1].cells);
          std::cout << std::setw(10) << 2.*log(rows[i-1].eta/rows[i].eta)/log(ratio);
        }
      std::cout << std::endl;
    }

  std::cout << std::defaultfloat << std::setprecision(6);
}

int main()
{
  // The strength of the coefficient is decided here, at run time
  const double contrast = 5.;

  Diffusion a(contrast);
  RightHandSide f(a);

  GenericEstimator Estimator;
  Estimator.AddVolumeTerm(WeightedVolumeResidual(a));
  Estimator.AddEdgeTerm(WeightedEdgeJump(a), INTERIOR_EDGES);

  GenericEstimator Energy;
  Energy.AddVolumeTerm(EnergyError(a));

  std::vector<Row> rows;

  UnitSquareMesh mesh(2);

  for(int level=0; level<6; ++level)
    {
      LagrangeElement element(1);
      FESpace Space(mesh, element);

      // The very same functors the estimator holds drive the discretization
      BilinearForm A(Space, Space);
      A.AddDiffusionTerm(a);
      A.Assemble();

      LinearForm F(Space);
      F.AddVolumeForce(f);
      F.Assemble();

      FEFunction Sol(Space);
      Sol.CreateFunction(A.SystemMatrix().Solve(F.LoadVector()));

      Row row;
      row.cells = mesh.NrCells();
      row.eta   = sqrt(Sum(Estimator.Assemble(Sol)));
      row.error = sqrt(Sum(Energy.Assemble(Sol)));
      rows.push_back(row);

      mesh.RefineUniform();
      mesh.RefineUniform();
    }

  std::cout << "Diffusionsproblem  -div(a grad u) = f,  a = 1 + " << contrast << "*x*y\n";
  std::cout << "Koeffizient und rechte Seite sind Funktoren auf Coordinate.\n";

  PrintTable("Uniforme Verfeinerung", rows);

  std::cout << "\nDiffusionTest was successful.\n";
  return 0;
}
