#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

#include "fem/BlockSystem.h"
#include "fem/DirichletValues.h"
#include "fem/LagrangeElement.h"
#include "fem/VtkOutput.h"
#include "mesh/Mesh.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// Flow through a channel with a cylinder, read from a mesh of gmsh. The geometry is the one
// of the DFG benchmark: a channel of 2.2 x 0.41 with a cylinder of radius 0.05 around
// (0.2, 0.2), see meshes/cylinder.geo.
//
//   du/dt - nu Laplace(u) + grad(p) = 0,  div(u) = 0
//
// with Taylor-Hood P2/P1 in space and the implicit Euler method in time. The fluid starts at
// rest and the parabolic inflow profile is switched on at t = 0. Walls and cylinder are no
// slip, at the outflow nothing is prescribed. That is the natural condition of the weak
// form, and it also determines the pressure, so no normalization is needed there.
//
// Every second step is written to cylinder_0001.vtk, cylinder_0002.vtk, ...

const double Length = 2.2;
const double Height = 0.41;
const double MaxInflow = 0.3;

/// Parabolic profile of the inflow, its mean value is 2/3 of MaxInflow
double Inflow(const Coordinate& p)
{
  return 4.*MaxInflow*p.y*(Height - p.y)/(Height*Height);
}

/// Everything except the outflow on the right, so walls, cylinder and inflow
bool NoOutflow(const Coordinate& p)
{
  return p.x < Length - 1.e-8;
}

bool InflowBoundary(const Coordinate& p)
{
  return p.x < 1.e-8;
}

bool Nowhere(const Coordinate&)
{
  return false;
}

/// (u, v)/tau + nu (grad u, grad v), the part of a time step acting on the new velocity
struct ViscousStep
{
  double nu, tau;

  double operator()(const PointValues& u, const PointValues& v) const
  {
    return u.value*v.value/tau + nu*dot(u.gradient, v.gradient);
  }
};

// -(du/dx, q) and -(du/dy, q). The transposed blocks give the pressure terms -(p, div v).
double DivergenceX(const PointValues& u, const PointValues& q) { return -u.gradient.x * q.value; }
double DivergenceY(const PointValues& u, const PointValues& q) { return -u.gradient.y * q.value; }

/// (u_old, v)/tau, with the velocity of the previous step as an FE function
struct OldVelocity
{
  const FEFunction& Uold;
  double tau;

  double operator()(const QuadPoint& p, const PointValues& v) const
  {
    return Uold.Value(p)/tau * v.value;
  }
};

/**
 * The implicit Euler method for the Stokes equations with the viscosity nu and the step size
 * tau. The velocity components need one bilinear form each, because they carry different
 * prescribed values.
 */
class ImplicitEuler
{
public:
  ImplicitEuler(Mesh& mesh, double nu, double tau)
    : P2(2), P1(1), V(mesh, P2, NoOutflow), Q(mesh, P1, Nowhere),
      Ux(V), Uy(V), P(Q), gx(V), gy(V),
      Ax(V, V), Ay(V, V), Bx(V, Q), By(V, Q), Fx(V), Fy(V),
      S({V, V, Q}), tau(tau)
  {
    gx.Set(Inflow, InflowBoundary);      // gy stays zero everywhere

    Ax.AddVolumeTerm(ViscousStep{nu, tau});
    Ay.AddVolumeTerm(ViscousStep{nu, tau});

    Bx.AddVolumeTerm(DivergenceX);
    By.AddVolumeTerm(DivergenceY);

    Fx.AddVolumeTerm(OldVelocity{Ux, tau});
    Fy.AddVolumeTerm(OldVelocity{Uy, tau});

    S.AddBlock(0, 0, Ax);
    S.AddBlock(1, 1, Ay);
    S.AddBlock(2, 0, Bx);
    S.AddBlock(2, 1, By);
    S.AddTransposedBlock(0, 2, Bx);
    S.AddTransposedBlock(1, 2, By);
    S.AddRhs(0, Fx);
    S.AddRhs(1, Fy);
    S.SetDirichletValues(0, gx);
    S.SetDirichletValues(1, gy);

    S.AssembleMatrix();
  }

  void Step()
  {
    S.AssembleRhs();

    // Many right hand sides with the same matrix, so the iterative refinement of UMFPACK is
    // not worth its three-fold cost here
    const Vector X = S.Solve(false);

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

  size_t NrDof() const { return S.NrDof(); }

  LagrangeElement P2, P1;
  FESpace V, Q;

  FEFunction Ux, Uy, P;

private:
  DirichletValues gx, gy;

  BilinearForm Ax, Ay, Bx, By;
  LinearForm Fx, Fy;
  BlockSystem S;

  const double tau;
};

int main()
{
  const double nu = 0.01;
  const double T = 2.;
  const int steps = 100;

  Mesh mesh("meshes/cylinder.msh");

  if(mesh.NrCells() == 0)
    {
      std::cerr << "ERROR: the mesh could not be read. Generate it with "
                << "gmsh -2 meshes/cylinder.geo -o meshes/cylinder.msh\n";
      return 1;
    }

  ImplicitEuler Euler(mesh, nu, T/steps);

  std::cout << "Channel with a cylinder, " << mesh.NrCells() << " cells, "
            << Euler.NrDof() << " unknowns, " << steps << " steps up to t = " << T
            << std::endl;

  VtkOutput out(mesh);
  out.AddVector("u", Euler.Ux, Euler.Uy);
  out.AddScalar("p", Euler.P);

  for(int n=1; n<=steps; ++n)
    {
      Euler.Step();

      if(n % 2 == 0)
        {
          std::ostringstream name;
          name << "cylinder_" << std::setw(4) << std::setfill('0') << n/2 << ".vtk";

          out.Write(name.str());
        }

      if(n % 20 == 0)
        std::cout << std::fixed << std::setprecision(2) << "  t = " << n*T/steps
                  << "   max speed " << std::setprecision(4) << Euler.MaxSpeed() << std::endl;
    }

  const double speed = Euler.MaxSpeed();

  // The inflow reaches 0.3, the fluid accelerates beside the cylinder, so the maximum has to
  // settle somewhat above that
  if(!(speed > 0.3 && speed < 1.))
    {
      std::cerr << "ERROR: the flow does not develop as expected, maximum speed " << speed
                << ".\n";
      return 1;
    }

  std::cout << "\nInstationaryStokesTest was successful, " << steps/2
            << " frames written to cylinder_*.vtk.\n";
  return 0;
}
