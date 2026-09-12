#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "fem/BlockSystem.h"
#include "fem/LagrangeElement.h"
#include "fem/ErrorNorm.h"
#include "fem/VtkOutput.h"
#include "mesh/UnitSquareMesh.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// Stokes equations  -Laplace(u) + grad(p) = f,  div(u) = 0  on the unit square with
// u = 0 on the boundary, discretized with the Taylor-Hood pair P2/P1. The velocity is
// the curl of the stream function x^2 (1-x)^2 y^2 (1-y)^2, the pressure has mean value
// zero, which the discrete pressure gets from a Lagrange multiplier.

// g(t) = t^2 (1-t)^2 and its derivatives
double g0(double t) { return t*t*(1.-t)*(1.-t); }
double g1(double t) { return 2.*t*(1.-t)*(1.-2.*t); }
double g2(double t) { return 2.*(1. - 6.*t + 6.*t*t); }
double g3(double t) { return 12.*(2.*t - 1.); }

double ux(const Coordinate& p) { return  g0(p.x)*g1(p.y); }
double uy(const Coordinate& p) { return -g1(p.x)*g0(p.y); }

Vector2D grad_ux(const Coordinate& p) { return Vector2D( g1(p.x)*g1(p.y),  g0(p.x)*g2(p.y)); }
Vector2D grad_uy(const Coordinate& p) { return Vector2D(-g2(p.x)*g0(p.y), -g1(p.x)*g1(p.y)); }

double pressure(const Coordinate& p)
{
  return cos(M_PI*p.x)*cos(M_PI*p.y);
}

double fx(const Coordinate& p)
{
  const double laplace = g2(p.x)*g1(p.y) + g0(p.x)*g3(p.y);
  return -laplace - M_PI*sin(M_PI*p.x)*cos(M_PI*p.y);
}

double fy(const Coordinate& p)
{
  const double laplace = -(g3(p.x)*g0(p.y) + g1(p.x)*g2(p.y));
  return -laplace - M_PI*cos(M_PI*p.x)*sin(M_PI*p.y);
}

// -(du/dx, q) and -(du/dy, q). The transposed blocks give the pressure terms -(p, div v).
double DivergenceX(const PointValues& u, const PointValues& q) { return -u.gradient.x * q.value; }
double DivergenceY(const PointValues& u, const PointValues& q) { return -u.gradient.y * q.value; }

bool Nowhere(const Coordinate&)
{
  return false;
}

double Eoc(const std::vector<double>& e, size_t i)
{
  return log2(e[i-1]/e[i]);
}

int main()
{
  const int levels = 5;

  Mesh mesh = UnitSquareMesh(3);

  std::vector<double> h1_u, l2_u, l2_p;

  std::cout << std::setw(8) << "Cells" << std::setw(8) << "DOFs"
            << std::setw(12) << "u H1" << std::setw(8) << "eoc"
            << std::setw(12) << "u L2" << std::setw(8) << "eoc"
            << std::setw(12) << "p L2" << std::setw(8) << "eoc" << std::endl;

  for(int level=0; level<levels; ++level)
    {
      LagrangeElement P2(2), P1(1);

      FESpace V(mesh, P2);
      FESpace Q(mesh, P1, Nowhere);

      BilinearForm A(V, V);
      A.AddLaplaceTerm();

      // -(div u, q), split into the two directions
      BilinearForm Bx(V, Q), By(V, Q);
      Bx.AddVolumeTerm(DivergenceX);
      By.AddVolumeTerm(DivergenceY);

      LinearForm Fx(V), Fy(V);
      Fx.AddVolumeForce(fx);
      Fy.AddVolumeForce(fy);

      BlockSystem S({V, V, Q});
      S.AddBlock(0, 0, A);
      S.AddBlock(1, 1, A);
      S.AddBlock(2, 0, Bx);
      S.AddBlock(2, 1, By);
      S.AddTransposedBlock(0, 2, Bx);
      S.AddTransposedBlock(1, 2, By);
      S.AddRhs(0, Fx);
      S.AddRhs(1, Fy);
      S.FixDof(2);
      S.Assemble();

      const Vector X = S.Solve();

      FEFunction Ux = S.Extract(0, X);
      FEFunction Uy = S.Extract(1, X);
      FEFunction P = S.Extract(2, X);

      // The fixed DOF leaves the pressure with an arbitrary constant, the exact one has
      // mean value zero
      P.SubtractMean();

      ErrorNorm Ex(ux, grad_ux), Ey(uy, grad_uy), Ep(pressure);

      h1_u.push_back(hypot(Ex.Compute(Ux, H1_SEMI), Ey.Compute(Uy, H1_SEMI)));
      l2_u.push_back(hypot(Ex.Compute(Ux, L2), Ey.Compute(Uy, L2)));
      l2_p.push_back(Ep.Compute(P, L2));

      std::cout << std::setw(8) << mesh.NrCells() << std::setw(8) << S.NrDof();
      const std::vector<double>* columns[3] = {&h1_u, &l2_u, &l2_p};
      for(int c=0; c<3; ++c)
        {
          std::cout << std::scientific << std::setprecision(3) << std::setw(12)
                    << columns[c]->back() << std::fixed << std::setw(8);
          if(level > 0)
            std::cout << Eoc(*columns[c], level);
          else
            std::cout << "";
        }
      std::cout << std::endl;

      if(level+1 == levels)
        {
          VtkOutput out(mesh);
          out.AddVector("u", Ux, Uy);
          out.AddScalar("p", P);
          out.Write("stokes.vtk");
        }

      if(level+1 < levels)
        {
          mesh.RefineUniform();
          mesh.RefineUniform();
        }
    }

  bool ok = true;

  if(std::fabs(Eoc(h1_u, levels-1) - 2.) > 0.15)
    {
      std::cerr << "ERROR: the H1 error of the velocity does not converge with order 2.\n";
      ok = false;
    }

  if(std::fabs(Eoc(l2_u, levels-1) - 3.) > 0.15)
    {
      std::cerr << "ERROR: the L2 error of the velocity does not converge with order 3.\n";
      ok = false;
    }

  if(std::fabs(Eoc(l2_p, levels-1) - 2.) > 0.15)
    {
      std::cerr << "ERROR: the L2 error of the pressure does not converge with order 2.\n";
      ok = false;
    }

  if(!ok)
    return 1;

  std::cout << "\nStokesTest was successful.\n";
  return 0;
}
