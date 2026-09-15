#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "fem/BlockSystem.h"
#include "fem/RaviartThomasElement.h"
#include "fem/DGElement.h"
#include "fem/ErrorNorm.h"
#include "fem/VtkOutput.h"
#include "mesh/UnitSquareMesh.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// Poisson in its mixed form on the unit square, with sigma = -grad(u) as a second unknown:
//
//   (sigma, tau) - (u, div tau) = -<g, tau.n>   for all tau in RT_0
//   -(div sigma, v)             = -(f, v)       for all v in DG_0
//
// The flux lives in the lowest order Raviart-Thomas space, the solution in the piecewise
// constants. That pair is inf-sup stable because div(RT_0) = DG_0 holds exactly. The
// Dirichlet data of u is natural here, it enters the right hand side over the boundary, so
// neither space carries an essential condition. Both unknowns converge with order 1.

double exact(const Coordinate& p)
{
  return sin(M_PI*p.x)*cos(M_PI*p.y);
}

/// -Laplace of the exact solution
double force(const Coordinate& p)
{
  return 2.*M_PI*M_PI*exact(p);
}

/// The two components of the exact flux sigma = -grad(u)
double flux_x(const Coordinate& p)
{
  return -M_PI*cos(M_PI*p.x)*cos(M_PI*p.y);
}

double flux_y(const Coordinate& p)
{
  return M_PI*sin(M_PI*p.x)*sin(M_PI*p.y);
}

bool Nowhere(const Coordinate&)
{
  return false;
}

int main()
{
  const int levels = 5;

  Mesh mesh = UnitSquareMesh(3);

  std::vector<double> u_errors, sigma_errors;
  bool ok = true;

  std::cout << std::setw(8) << "Cells" << std::setw(10) << "DOFs"
            << std::setw(12) << "u L2" << std::setw(8) << "eoc"
            << std::setw(12) << "sigma L2" << std::setw(8) << "eoc" << std::endl;

  for(int level=0; level<levels; ++level)
    {
      RaviartThomasElement rt;
      DGElement dg(0);

      FESpace V(mesh, rt, Nowhere), Q(mesh, dg, Nowhere);

      if(V.NrDof() != mesh.NrEdges() || Q.NrDof() != mesh.NrCells())
        {
          std::cerr << "ERROR: the spaces have " << V.NrDof() << " and " << Q.NrDof()
                    << " DOFs, the mesh has " << mesh.NrEdges() << " edges and "
                    << mesh.NrCells() << " cells.\n";
          ok = false;
        }

      BilinearForm A(V, V), B(V, Q);

      A.AddVolumeTerm([](const VectorValues& sigma, const VectorValues& tau)
                      { return dot(sigma.value, tau.value); });

      // The transposed block gives -(u, div tau)
      B.AddVolumeTerm([](const VectorValues& sigma, const PointValues& v)
                      { return -sigma.divergence * v.value; });

      // The Dirichlet data of u is natural in the mixed form
      LinearForm F(V);
      F.AddBoundaryTerm([](const QuadPoint& p, const EdgeGeometry& edge,
                           const VectorValues& tau)
                        { return -exact(p.x) * dot(tau.value, edge.normal); });

      LinearForm G(Q);
      G.AddVolumeTerm([](const QuadPoint& p, const PointValues& v)
                      { return -force(p.x) * v.value; });

      BlockSystem S({V, Q});
      S.AddBlock(0, 0, A);
      S.AddBlock(1, 0, B);
      S.AddTransposedBlock(0, 1, B);
      S.AddRhs(0, F);
      S.AddRhs(1, G);

      S.AssembleMatrix();

      const Vector X = S.Solve(S.AssembleRhs());

      FEFunction Sigma = S.Extract(0, X);
      FEFunction U = S.Extract(1, X);

      ErrorNorm Eu(exact), Ex(flux_x), Ey(flux_y);

      u_errors.push_back(Eu.Compute(U, L2));
      sigma_errors.push_back(std::hypot(Ex.Compute(Sigma, L2, 0), Ey.Compute(Sigma, L2, 1)));

      std::cout << std::setw(8) << mesh.NrCells() << std::setw(10) << S.NrDof()
                << std::scientific << std::setprecision(3) << std::setw(12)
                << u_errors.back();
      if(level > 0)
        std::cout << std::fixed << std::setw(8)
                  << log2(u_errors[level-1]/u_errors[level]);
      else
        std::cout << std::setw(8) << "";
      std::cout << std::scientific << std::setw(12) << sigma_errors.back();
      if(level > 0)
        std::cout << std::fixed << std::setw(8)
                  << log2(sigma_errors[level-1]/sigma_errors[level]);
      std::cout << std::endl;

      if(level+1 == levels)
        {
          VtkOutput out(mesh);
          out.AddVector("sigma", Sigma);
          out.AddScalar("u", U);
          out.Write("mixed_poisson.vtk");
        }

      if(level+1 < levels)
        {
          mesh.RefineUniform();
          mesh.RefineUniform();
        }
    }

  const double u_eoc = log2(u_errors[levels-2]/u_errors[levels-1]);
  const double sigma_eoc = log2(sigma_errors[levels-2]/sigma_errors[levels-1]);

  if(std::fabs(u_eoc - 1.) > 0.15)
    {
      std::cerr << "ERROR: the error of u does not converge with order 1.\n";
      ok = false;
    }

  if(std::fabs(sigma_eoc - 1.) > 0.15)
    {
      std::cerr << "ERROR: the error of sigma does not converge with order 1.\n";
      ok = false;
    }

  if(!ok)
    return 1;

  std::cout << "\nMixedPoissonTest was successful.\n";
  return 0;
}
