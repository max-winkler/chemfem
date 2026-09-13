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

/// Largest speed in the degrees of freedom of a vector valued function
double MaxSpeed(const FEFunction& u)
{
  double speed = 0.;
  for(size_t k=0; k+1<u.GetFESpace().NrDof(); k+=2)
    speed = std::max(speed, hypot(u[k], u[k+1]));
  return speed;
}

int main()
{
  const double nu = 0.01;
  const double T = 2.;
  const int steps = 100;
  const double tau = T/steps;

  Mesh mesh("meshes/cylinder.msh");

  if(mesh.NrCells() == 0)
    {
      std::cerr << "ERROR: the mesh could not be read. Generate it with "
                << "gmsh -2 meshes/cylinder.geo -o meshes/cylinder.msh\n";
      return 1;
    }

  LagrangeElement P2(2), P1(1);
  ProductElement Velocity(P2, 2);

  FESpace V(mesh, Velocity, NoOutflow), Q(mesh, P1, Nowhere);

  DirichletValues g(V);
  g.Set(Inflow, InflowBoundary);

  BilinearForm A(V, V), B(V, Q), M(V, V);

  A.AddVolumeTerm([nu, tau](const VectorValues& u, const VectorValues& v)
                  { return dot(u.value, v.value)/tau + nu*ddot(u.gradient, v.gradient); });

  // The transposed block gives the pressure term -(p, div v)
  B.AddVolumeTerm([](const VectorValues& u, const PointValues& q)
                  { return -u.divergence * q.value; });

  M.AddVolumeTerm([](const VectorValues& u, const VectorValues& v)
                  { return dot(u.value, v.value); });

  BlockSystem S({V, Q});
  S.AddBlock(0, 0, A);
  S.AddBlock(1, 0, B);
  S.AddTransposedBlock(0, 1, B);
  S.SetDirichletValues(0, g);
  S.AssembleMatrix();

  // The old velocity enters the right hand side as (1/tau) M u_old, its Dirichlet columns
  // hold the inflow profile
  M.SetDirichletValues(g);
  M.Assemble();

  // Constant in time, only the mass term is added in each step
  const Vector Rhs0 = S.AssembleRhs();

  FEFunction U(V), P(Q);

  // Velocity of the previous step, free DOFs only
  Vector UFree(V.NrFreeDof());

  std::cout << "Channel with a cylinder, " << mesh.NrCells() << " cells, "
            << S.NrDof() << " unknowns, " << steps << " steps up to t = " << T
            << std::endl;

  VtkOutput out(mesh);
  out.AddVector("u", U);
  out.AddScalar("p", P);

  for(int n=1; n<=steps; ++n)
    {
      Vector Rhs = Rhs0;
      S.AddToRhs(Rhs, 0, (1./tau) * (M.SystemMatrix()*UFree + M.DirichletRhs()));

      const Vector X = S.Solve(Rhs);

      UFree = S.FreeDof(0, X);

      U = S.Extract(0, X);
      P = S.Extract(1, X);

      if(n % 2 == 0)
        {
          std::ostringstream name;
          name << "cylinder_" << std::setw(4) << std::setfill('0') << n/2 << ".vtk";

          out.Write(name.str());
        }

      if(n % 20 == 0)
        std::cout << std::fixed << std::setprecision(2) << "  t = " << n*T/steps
                  << "   max speed " << std::setprecision(4) << MaxSpeed(U) << std::endl;
    }

  const double speed = MaxSpeed(U);

  // The inflow reaches 0.3 and the fluid accelerates beside the cylinder
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
