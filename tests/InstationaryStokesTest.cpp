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

// Flow through a channel with a cylinder, read from a mesh of gmsh. The geometry is the one
// of the DFG benchmark: a channel of 2.2 x 0.41 with a cylinder of radius 0.05 around
// (0.2, 0.2), see meshes/cylinder.geo.
//
//   du/dt - nu Laplace(u) + grad(p) = 0,  div(u) = 0
//
// with Taylor-Hood P2/P1 in space and the implicit Euler method in time. The velocity lives
// in one vector valued space. The fluid starts at rest and the parabolic inflow profile is
// switched on at t = 0. Walls and cylinder are no slip, at the outflow nothing is
// prescribed. That is the natural condition of the weak form, and it also determines the
// pressure, so no normalization is needed there.
//
// Every second step is written to cylinder_0001.vtk, cylinder_0002.vtk, ...

const double Length = 2.2;
const double Height = 0.41;
const double MaxInflow = 0.3;

/// Parabolic profile of the inflow, its mean value is 2/3 of MaxInflow
Vector2D Inflow(const Coordinate& p)
{
  return Vector2D(4.*MaxInflow*p.y*(Height - p.y)/(Height*Height), 0.);
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

  double operator()(const VectorValues& u, const VectorValues& v) const
  {
    return dot(u.value, v.value)/tau + nu*ddot(u.gradient, v.gradient);
  }
};

/// -(div u, q). The transposed block gives the pressure term -(p, div v).
double Divergence(const VectorValues& u, const PointValues& q)
{
  return -u.divergence * q.value;
}

/// (u_old, v)/tau, with the velocity of the previous step as an FE function
struct OldVelocity
{
  const FEFunction& Uold;
  double tau;

  double operator()(const QuadPoint& p, const VectorValues& v) const
  {
    return dot(Uold.VectorValue(p), v.value)/tau;
  }
};

/**
 * The implicit Euler method for the Stokes equations with the viscosity nu and the step size
 * tau.
 */
class ImplicitEuler
{
public:
  ImplicitEuler(Mesh& mesh, double nu, double tau)
    : P2(2), P1(1), Velocity(P2, 2), V(mesh, Velocity, NoOutflow), Q(mesh, P1, Nowhere),
      U(V), P(Q), g(V), A(V, V), B(V, Q), F(V), S({V, Q}), tau(tau)
  {
    g.Set(Inflow, InflowBoundary);

    A.AddVolumeTerm(ViscousStep{nu, tau});
    B.AddVolumeTerm(Divergence);
    F.AddVolumeTerm(OldVelocity{U, tau});

    S.AddBlock(0, 0, A);
    S.AddBlock(1, 0, B);
    S.AddTransposedBlock(0, 1, B);
    S.AddRhs(0, F);
    S.SetDirichletValues(0, g);

    S.AssembleMatrix();
  }

  void Step()
  {
    S.AssembleRhs();

    // Many right hand sides with the same matrix, so the iterative refinement of UMFPACK is
    // not worth its three-fold cost here
    const Vector X = S.Solve(false);

    U = S.Extract(0, X);
    P = S.Extract(1, X);
  }

  /// Largest speed in the degrees of freedom
  double MaxSpeed() const
  {
    double speed = 0.;
    for(size_t k=0; k+1<V.NrDof(); k+=2)
      speed = std::max(speed, hypot(U[k], U[k+1]));
    return speed;
  }

  size_t NrDof() const { return S.NrDof(); }

  LagrangeElement P2, P1;
  ProductElement Velocity;
  FESpace V, Q;

  FEFunction U, P;

private:
  DirichletValues g;

  BilinearForm A, B;
  LinearForm F;
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
  out.AddVector("u", Euler.U);
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
