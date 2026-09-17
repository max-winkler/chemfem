#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

#include "fem/BlockSystem.h"
#include "fem/DirichletValues.h"
#include "fem/ProductElement.h"
#include "fem/LagrangeElement.h"
#include "fem/VtkOutput.h"
#include "mesh/Mesh.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// Flow through a channel with a cylinder, the geometry of the DFG benchmark, see
// meshes/cylinder.geo. In contrast to InstationaryStokesTest the convection is kept:
//
//   du/dt + (u.grad)u - nu Laplace(u) + grad(p) = 0,  div(u) = 0
//
// with Taylor-Hood P2/P1 in space and the implicit Euler method in time. Every time step is
// a nonlinear problem, solved with Newton's method. The increment carries homogeneous
// boundary values, which is why the Jacobian never gets SetDirichletValues: without
// prescribed values its lifting is zero. The iterate itself keeps the inflow profile,
// because CreateFunction writes the prescribed values into the Dirichlet DOFs of every
// iterate.
//
// With nu = 0.001 and a mean inflow of 1.0 over the cylinder diameter 0.1 the Reynolds
// number is 100, where the flow sheds vortices. The test checks that Newton converges in
// every step and that the solution stays bounded, and it writes the frames to
// navier_stokes_*.vtk. The drag and lift coefficients need a boundary functional, which the
// library does not have yet.
//
// The integrands capture the iterate by reference, so Velocity has to be declared before
// the forms that use it and must not be replaced, only refreshed.

const double Height = 0.41;
const double MaxInflow = 1.5;

const double nu = 0.001;
const double T = 1.0;
const int steps = 200;

const double NewtonTolerance = 1.e-10;
const int MaxNewton = 12;

/// Parabolic profile of the inflow, its mean value is 2/3 of MaxInflow
Vector2D Inflow(const Coordinate& p)
{
  return Vector2D(4.*MaxInflow*p.y*(Height - p.y)/(Height*Height), 0.);
}

/// Everything except the outflow on the right, so walls, cylinder and inflow
bool NoOutflow(const Coordinate& p)
{
  return p.x < 2.2 - 1.e-8;
}

bool InflowBoundary(const Coordinate& p)
{
  return p.x < 1.e-8;
}

/// (a.grad)w, the component i is the gradient of w_i times a
Vector2D Convect(const Matrix2D& grad_w, const Vector2D& a)
{
  return grad_w * a;
}

/// Largest speed in the degrees of freedom of a vector valued function
double MaxSpeed(const FEFunction& u)
{
  double speed = 0.;
  for(size_t k=0; k+1<u.GetFESpace().NrDof(); k+=2)
    speed = std::max(speed, std::hypot(u[k], u[k+1]));
  return speed;
}

int main()
{
  const double tau = T/steps;

  Mesh mesh("meshes/cylinder.msh");

  if(mesh.NrCells() == 0)
    {
      std::cerr << "ERROR: the mesh could not be read. Generate it with "
                << "gmsh -2 meshes/cylinder.geo -o meshes/cylinder.msh\n";
      return 1;
    }

  LagrangeElement P2(2), P1(1);
  ProductElement VelocityElement(P2, 2);

  FESpace V(mesh, VelocityElement, NoOutflow), Q(mesh, P1);

  DirichletValues g(V);
  g.Set(Inflow, InflowBoundary);

  // The iterate of the Newton method. The integrands below hold a reference to it, so it is
  // declared before them and refreshed in place.
  FEFunction Velocity(V), Pressure(Q);

  BilinearForm A(V, V), B(V, Q), M(V, V);

  // The linear part of the Jacobian, as in the Stokes case
  A.AddVolumeTerm([tau](const VectorValues& u, const VectorValues& v)
                  { return dot(u.value, v.value)/tau + nu*ddot(u.gradient, v.gradient); });

  // The convection, linearized around the current iterate
  A.AddVolumeTerm([&Velocity](const QuadPoint& p, const VectorValues& u,
                              const VectorValues& v)
                  {
                    const VectorValues w = Velocity.EvaluateVector(p);

                    return dot(Convect(u.gradient, w.value), v.value)
                         + dot(Convect(w.gradient, u.value), v.value);
                  });

  // The transposed block gives the pressure term -(p, div v)
  B.AddVolumeTerm([](const VectorValues& u, const PointValues& q)
                  { return -u.divergence * q.value; });

  M.AddVolumeTerm([](const VectorValues& u, const VectorValues& v)
                  { return dot(u.value, v.value); });

  // The residual of the momentum equation, without the terms that the matrix and the mass
  // term already cover. The sign is the one of the right hand side.
  LinearForm F(V);
  F.AddVolumeTerm([&Velocity, &Pressure](const QuadPoint& p, const VectorValues& v)
                  {
                    const VectorValues w = Velocity.EvaluateVector(p);
                    const double q = Pressure.Value(p);

                    return -dot(Convect(w.gradient, w.value), v.value)
                         - nu*ddot(w.gradient, v.gradient)
                         + q*v.divergence;
                  });

  // The residual of the continuity equation
  LinearForm G(Q);
  G.AddVolumeTerm([&Velocity](const QuadPoint& p, const PointValues& q)
                  { return Velocity.EvaluateVector(p).divergence * q.value; });

  BlockSystem S({V, Q});
  S.AddBlock(0, 0, A);
  S.AddBlock(1, 0, B);
  S.AddTransposedBlock(0, 1, B);
  S.AddRhs(0, F);
  S.AddRhs(1, G);

  // No SetDirichletValues: the increment of the Newton method is homogeneous

  // The mass matrix is the same in every step. It carries no prescribed values either, the
  // Dirichlet part of the difference of two iterates vanishes.
  M.Assemble();

  VtkOutput out(mesh);
  out.AddVector("u", Velocity);
  out.AddScalar("p", Pressure);

  // The fluid starts at rest, the inflow is switched on at t = 0. CreateFunction puts the
  // prescribed values into the Dirichlet DOFs, so the initial state satisfies them exactly.
  Vector VelocityDof(V.NrFreeDof()), PressureDof(Q.NrFreeDof());

  Velocity.CreateFunction(VelocityDof, g);
  Pressure.CreateFunction(PressureDof);

  std::cout << "Navier-Stokes around a cylinder, " << mesh.NrCells() << " cells, "
            << S.NrDof() << " unknowns, Re = 100\n"
            << std::setw(8) << "step" << std::setw(8) << "Newton"
            << std::setw(14) << "last update" << std::setw(14) << "max speed" << std::endl;

  bool ok = true;

  for(int n=1; n<=steps; ++n)
    {
      const Vector Old = VelocityDof;

      // The time derivative enters through the mass matrix on the DOF vectors, which is
      // cheaper than integrating it again. It changes within the Newton loop, because the
      // iterate does.
      int iteration = 0;
      double update = 0.;

      for(; iteration < MaxNewton; ++iteration)
        {
          S.AssembleMatrix();

          Vector Rhs = S.AssembleRhs();
          S.AddToRhs(Rhs, 0, (1./tau) * (M.SystemMatrix() * (Old - VelocityDof)));

          const Vector X = S.Solve(Rhs);

          const Vector dVelocity = S.FreeDof(0, X);
          const Vector dPressure = S.FreeDof(1, X);

          update = dVelocity.Norm();

          // The updates of the first steps show that Newton converges quadratically, which
          // is what tells a correct Jacobian from a merely convergent one
          if(n <= 2)
            std::cout << std::setw(8) << n << std::setw(8) << iteration+1
                      << std::scientific << std::setprecision(3) << std::setw(14) << update
                      << std::endl;

          VelocityDof += dVelocity;
          PressureDof += dPressure;

          Velocity.CreateFunction(VelocityDof, g);
          Pressure.CreateFunction(PressureDof);

          if(update < NewtonTolerance)
            break;
        }

      const double speed = MaxSpeed(Velocity);

      if(iteration == MaxNewton)
        {
          std::cerr << "ERROR: Newton did not converge in step " << n << ", the last update "
                    << "was " << update << ".\n";
          ok = false;
          break;
        }

      if(!(speed < 10.*MaxInflow))
        {
          std::cerr << "ERROR: the velocity grew to " << speed << " in step " << n
                    << ", the solution is not bounded.\n";
          ok = false;
          break;
        }

      if(n % 10 == 0)
        std::cout << std::setw(8) << n << std::setw(8) << iteration+1
                  << std::scientific << std::setprecision(3) << std::setw(14) << update
                  << std::fixed << std::setprecision(4) << std::setw(14) << speed
                  << std::endl;

      if(n % 4 == 0)
        {
          std::ostringstream name;
          name << "navier_stokes_" << std::setw(4) << std::setfill('0') << n/4 << ".vtk";
          out.Write(name.str());
        }
    }

  if(!ok)
    return 1;

  std::cout << "\nNavierStokesTest was successful, " << steps/4
            << " frames written to navier_stokes_*.vtk.\n";
  return 0;
}
