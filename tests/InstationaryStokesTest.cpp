#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <memory>

#include "fem/BlockSystem.h"
#include "fem/LagrangeElement.h"
#include "fem/VtkOutput.h"
#include "mesh/UnitSquareMesh.h"
#include "linalg/DirectSolver.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// Stirred fluid in a closed box, described by the instationary Stokes equations
//
//   du/dt - nu Laplace(u) + grad(p) = f,  div(u) = 0,  u = 0 on the walls,
//
// with Taylor-Hood P2/P1 in space and the implicit Euler method in time. Each step solves
//
//   (u, v)/tau + nu (grad u, grad v) - (p, div v) = (f, v) + (u_old, v)/tau
//                                     -(div u, q) = 0
//
// whose matrix does not change, so it is factorized only once.
//
// The fluid starts at rest and is driven by two stirrers, forces with the profile of a
// Gaussian vortex, which circle around the center of the box. The low viscosity lets the
// fluid keep its momentum, so the whole box is spun up over the revolutions.
//
// Every second step is written to stokes_flow_0001.vtk, stokes_flow_0002.vtk, ... Open
// the series in ParaView and show u e.g. with the Surface LIC representation.

/**
 * Force of the two stirrers. The functor holds a reference to the time of the time loop,
 * so the copies the linear forms keep of it always see the current time.
 */
struct Stirring
{
  const double& t;
  int component;

  double operator()(const Coordinate& p) const
  {
    const double width = 0.005;

    double f = 0.;

    // Two stirrers on opposite sides of a circle of radius 0.25, one revolution per time unit
    for(int s=0; s<2; ++s)
      {
        const double angle = 2.*M_PI*t + s*M_PI;
        const double dx = p.x - (0.5 + 0.25*cos(angle));
        const double dy = p.y - (0.5 + 0.25*sin(angle));

        const double bump = exp(-(dx*dx + dy*dy)/width);

        // The force is the curl (d/dy, -d/dx) of the bump, a counterclockwise vortex
        f += component == 0 ? -2.*dy/width*bump : 2.*dx/width*bump;
      }

    return f;
  }
};

/// (u, v)/tau + nu (grad u, grad v), the part of a time step acting on the new velocity
struct ViscousStep
{
  double nu, tau;

  double operator()(const QuadPoint&, const CellGeometry&,
                    const PointValues& u, const PointValues& v) const
  {
    return u.value*v.value/tau + nu*dot(u.gradient, v.gradient);
  }
};

/// -(p, dv/dx_i), the trial function is the pressure
struct PressureTerm
{
  int direction;

  double operator()(const QuadPoint&, const CellGeometry&,
                    const PointValues& p, const PointValues& v) const
  {
    return -p.value * v.gradient[direction];
  }
};

/// -(du/dx_i, q), the test function belongs to the pressure
struct DivergenceTerm
{
  int direction;

  double operator()(const QuadPoint&, const CellGeometry&,
                    const PointValues& u, const PointValues& q) const
  {
    return -u.gradient[direction] * q.value;
  }
};

/// (u_old, v)/tau, with the velocity of the previous step as an FE function
struct OldVelocity
{
  const FEFunction& Uold;
  double tau;

  double operator()(const QuadPoint& p, const CellGeometry&, const PointValues& v) const
  {
    return Uold.Evaluate(p).value/tau * v.value;
  }
};

bool Nowhere(const Coordinate&)
{
  return false;
}

/**
 * The implicit Euler method for the Stokes equations with the viscosity nu and the step
 * size tau. The fluid starts at rest, each call of Step() advances the time t by tau.
 */
class ImplicitEuler
{
public:
  ImplicitEuler(Mesh& mesh, double nu, double tau, double& t,
                ScalarFunction fx, ScalarFunction fy)
    : P2(2), P1(1), V(mesh, P2), Q(mesh, P1, Nowhere), Ux(V), Uy(V), P(Q),
      A(V, V), BxT(Q, V), ByT(Q, V), Bx(V, Q), By(V, Q), Fx(V), Fy(V),
      S({&V, &V, &Q}), tau(tau), t(t)
  {
    A.AddVolumeTerm(ViscousStep{nu, tau});

    // -(p, div v) and -(div u, q)
    BxT.AddVolumeTerm(PressureTerm{0});
    ByT.AddVolumeTerm(PressureTerm{1});
    Bx.AddVolumeTerm(DivergenceTerm{0});
    By.AddVolumeTerm(DivergenceTerm{1});

    Fx.AddVolumeForce(fx);
    Fx.AddVolumeTerm(OldVelocity{Ux, tau});
    Fy.AddVolumeForce(fy);
    Fy.AddVolumeTerm(OldVelocity{Uy, tau});

    S.AddBlock(0, 0, A);
    S.AddBlock(1, 1, A);
    S.AddBlock(0, 2, BxT);
    S.AddBlock(1, 2, ByT);
    S.AddBlock(2, 0, Bx);
    S.AddBlock(2, 1, By);
    S.AddRhs(0, Fx);
    S.AddRhs(1, Fy);
    S.AddMeanValueConstraint(2);

    S.AssembleMatrix();
    LU = std::make_unique<DirectSolver>(S.SystemMatrix());
  }

  void Step()
  {
    t += tau;

    S.AssembleRhs();
    const Vector X = LU->Solve(S.Rhs());

    Ux = S.Extract(0, X);
    Uy = S.Extract(1, X);
    P = S.Extract(2, X);
  }

  /// Largest speed in the degrees of freedom
  double MaxSpeed() const
  {
    double speed = 0.;
    for(size_t k=0; k<V.NrDof(); ++k)
      speed = std::max(speed, hypot(Ux[k], Uy[k]));
    return speed;
  }

  LagrangeElement P2, P1;
  FESpace V, Q;

  FEFunction Ux, Uy, P;

private:
  BilinearForm A, BxT, ByT, Bx, By;
  LinearForm Fx, Fy;
  BlockSystem S;
  std::unique_ptr<DirectSolver> LU;

  const double tau;
  double& t;
};

int main()
{
  const double nu = 0.01;
  const double T = 2.;
  const int steps = 200;

  Mesh mesh = UnitSquareMesh(3);
  for(int i=0; i<8; ++i)
    mesh.RefineUniform();

  double t = 0.;
  ImplicitEuler Euler(mesh, nu, T/steps, t, Stirring{t, 0}, Stirring{t, 1});

  std::cout << "Stirred fluid on " << mesh.NrCells() << " cells, " << steps
            << " steps up to t = " << T << std::endl;

  for(int n=1; n<=steps; ++n)
    {
      Euler.Step();

      if(n % 2 == 0)
        {
          std::ostringstream name;
          name << "stokes_flow_" << std::setw(4) << std::setfill('0') << n/2 << ".vtk";

          VtkOutput out(mesh);
          out.AddVector("u", Euler.Ux, Euler.Uy);
          out.AddScalar("p", Euler.P);
          out.Write(name.str());
        }

      if(n % 20 == 0)
        std::cout << std::fixed << std::setprecision(2) << "  t = " << t
                  << "   max speed " << std::scientific << std::setprecision(3)
                  << Euler.MaxSpeed() << std::endl;
    }

  const double speed = Euler.MaxSpeed();

  if(!(speed > 1.e-3 && speed < 1.e3))
    {
      std::cerr << "ERROR: the fluid is at rest or the simulation blows up.\n";
      return 1;
    }

  std::cout << "\nInstationaryStokesTest was successful, " << steps/2
            << " frames written to stokes_flow_*.vtk.\n";
  return 0;
}
