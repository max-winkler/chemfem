#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "fem/LinearForm.h"
#include "fem/BilinearForm.h"
#include "fem/LagrangeElement.h"
#include "fem/FEFunction.h"
#include "fem/GenericEstimator.h"
#include "mesh/LShapeMesh.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// ---------------------------------------------------------------------------------
// -Laplace(u) = f  on the L shaped domain, u = 0 on the boundary.
//
// The right hand side is a pair of narrow Gaussian peaks, one positive and one
// negative, placed in the two arms of the L. Together with the reentrant corner at
// the origin, where the solution behaves like r^(2/3), the problem has three
// features at three different places and nothing worth resolving in between.
//
// The adaptive loop therefore produces a strongly graded mesh: dense at the peaks,
// dense at the corner, untouched at its initial size elsewhere. That is also what
// makes it a sharp regression test for the Newest-Vertex-Bisection closure in
// Mesh::Refine, since strongly graded meshes are what used to produce hanging nodes
// and corrupted edge-neighbor relations. Mesh::Check() runs after every step.
// ---------------------------------------------------------------------------------

static double Peak(const Coordinate& p, double x0, double y0, double width)
{
  const double dx = p.x - x0, dy = p.y - y0;

  return exp(-width*(dx*dx + dy*dy));
}

double f(const Coordinate& p)
{
  return 60.*Peak(p, 0.55, 0.55, 300.) - 40.*Peak(p, -0.6, 0.45, 150.);
}

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

/**
 * Doerfler marking: mark the smallest set of cells carrying the fraction theta of
 * the total indicator. Unlike the maximum strategy this refines enough cells per
 * step to reach the optimal rate.
 */
std::vector<bool> MarkCells(const Vector& Indicators, double theta)
{
  std::vector<std::pair<double, size_t> > sorted;
  sorted.reserve(Indicators.size());

  size_t c = 0;
  for(Vector::const_iterator it = Indicators.begin(); it != Indicators.end(); ++it, ++c)
    sorted.push_back(std::make_pair(*it, c));

  std::sort(sorted.begin(), sorted.end(),
            std::greater<std::pair<double, size_t> >());

  const double target = theta * Sum(Indicators);

  std::vector<bool> Marker(Indicators.size(), false);
  double collected = 0.;

  for(size_t i=0; i<sorted.size() && collected < target; ++i)
    {
      Marker[sorted[i].second] = true;
      collected += sorted[i].first;
    }

  return Marker;
}

/// Ratio of the largest to the smallest cell diameter
double Grading(const Mesh& mesh)
{
  double h_min = -1., h_max = 0.;

  for(size_t c=0; c<mesh.NrCells(); ++c)
    {
      const double h = mesh.GetCellInfo(c).Diam();

      if(h > h_max) h_max = h;
      if(h_min < 0. || h < h_min) h_min = h;
    }

  return h_max/h_min;
}

struct Row
{
  size_t cells, dofs;
  double eta;
};

/// Rate s in eta ~ N^(-s) between two consecutive rows
double Rate(const Row& a, const Row& b)
{
  return log(a.eta/b.eta) / log(double(b.dofs)/double(a.dofs));
}

void PrintTable(const std::string& caption, const std::vector<Row>& rows)
{
  std::cout << "\n" << caption << "\n" << std::string(43, '=') << "\n";
  std::cout << std::setw(9) << "Zellen" << std::setw(9) << "DOFs"
            << std::setw(14) << "eta" << std::setw(9) << "Rate" << std::endl;

  for(size_t i=0; i<rows.size(); ++i)
    {
      std::cout << std::setw(9) << rows[i].cells
                << std::setw(9) << rows[i].dofs
                << std::setw(14) << std::scientific << std::setprecision(4) << rows[i].eta;

      if(i > 0)
        std::cout << std::setw(9) << std::fixed << std::setprecision(3)
                  << Rate(rows[i-1], rows[i]);

      std::cout << std::endl;
    }

  std::cout << std::defaultfloat << std::setprecision(6);
}

int main()
{
  GenericEstimator Estimator = GenericEstimator::Residual(f);

  // ------------------------------------------------------------------ uniform
  std::vector<Row> uniform;

  {
    LShapeMesh mesh(4);

    for(int level=0; level<5; ++level)
      {
        LagrangeElement element(1);
        FESpace Space(mesh, element);

        FEFunction Sol = Solve(Space);

        Row row;
        row.cells = mesh.NrCells();
        row.dofs  = Space.NrFreeDof();
        row.eta   = sqrt(Sum(Estimator.Assemble(Sol)));
        uniform.push_back(row);

        if(level+1 < 5)
          {
            mesh.RefineUniform();
            mesh.RefineUniform();
          }
      }
  }

  // ----------------------------------------------------------------- adaptive
  std::vector<Row> adaptive;
  double grading = 1.;

  {
    LShapeMesh mesh(4);

    const int max_iter = 18;

    for(int level=0; level<max_iter; ++level)
      {
        LagrangeElement element(1);
        FESpace Space(mesh, element);

        FEFunction Sol = Solve(Space);

        Vector Indicators = Estimator.Assemble(Sol);

        Row row;
        row.cells = mesh.NrCells();
        row.dofs  = Space.NrFreeDof();
        row.eta   = sqrt(Sum(Indicators));
        adaptive.push_back(row);

        Sol.WriteVtk("adaptive_" + std::to_string(level) + ".vtk");

        if(level+1 == max_iter)
          {
            grading = Grading(mesh);
            break;
          }

        mesh.Refine(MarkCells(Indicators, 0.6));

        if(!mesh.Check())
          {
            std::cerr << "ERROR: the mesh is broken after refinement step "
                      << level+1 << ".\n";
            return 1;
          }
      }
  }

  PrintTable("Uniforme Verfeinerung", uniform);
  PrintTable("Adaptive Verfeinerung (Doerfler, theta = 0.6)", adaptive);

  // The adaptive row closest in size to the finest uniform one
  const Row& reference = uniform.back();

  size_t best = 0;
  for(size_t i=1; i<adaptive.size(); ++i)
    if(std::labs(long(adaptive[i].dofs) - long(reference.dofs))
       < std::labs(long(adaptive[best].dofs) - long(reference.dofs)))
      best = i;

  const double rate = Rate(adaptive[adaptive.size()-2], adaptive.back());
  const double gain = reference.eta / adaptive[best].eta;

  std::cout << "\n" << std::fixed << std::setprecision(3)
            << "Rate adaptiv                 : " << rate
            << "   (erwartet 1/2, optimal fuer P1)\n"
            << "eta uniform / eta adaptiv    : " << gain
            << "   bei " << reference.dofs << " gegen "
            << adaptive[best].dofs << " DOFs\n"
            << "Netzgradierung h_max / h_min : " << grading << std::endl;

  bool ok = true;

  if(std::fabs(rate - 0.5) > 0.1)
    {
      std::cerr << "ERROR: the adaptive rate is not the expected 1/2.\n";
      ok = false;
    }

  if(gain < 2.)
    {
      std::cerr << "ERROR: adaptive refinement does not pay off against uniform.\n";
      ok = false;
    }

  if(grading < 20.)
    {
      std::cerr << "ERROR: the mesh is not graded, the refinement spread out.\n";
      ok = false;
    }

  if(!ok) return 1;

  std::cout << "\nAdaptivityTest was successful.\n";
  return 0;
}
