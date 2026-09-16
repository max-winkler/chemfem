#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "Verbosity.h"
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
// -Laplace(u) = f  on the L shaped domain, u = 0 on the boundary, with two narrow
// Gaussian peaks in the arms of the L and the r^(2/3) corner at the origin.
//
// The strongly graded meshes this produces are a regression test for the closure in
// Mesh::Refine, which used to leave hanging nodes. From P2 on the test also covers
// the Laplacian in the volume residual.
// ---------------------------------------------------------------------------------

/// Polynomial degree of the Lagrange element
const int Degree = 3;

static_assert(Degree >= 1 && Degree <= 4,
              "LagrangeElement implements the degrees 1 to 4");

/// Optimal rate of the adaptive loop, eta ~ N^(-Degree/2)
const double OptimalRate = 0.5*Degree;

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

  Vector X(A.SystemMatrix().Solve(F.LoadVector(), LIN_SOLVER::UMFPACK));

  FEFunction Sol(Space);
  Sol.CreateFunction(X);

  return Sol;
}

/// Doerfler marking: the smallest set of cells carrying the fraction theta of the
/// total indicator
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

/// Rate s in eta ~ N^(-s) between two rows, N the number of degrees of freedom
double Rate(const Row& a, const Row& b)
{
  return log(a.eta/b.eta) / log(double(b.dofs)/double(a.dofs));
}

void PrintTable(const std::string& caption, const std::vector<Row>& rows)
{
  std::cout << "\n" << caption << "\n" << std::string(43, '=') << "\n";
  std::cout << std::setw(9) << "Cells" << std::setw(9) << "DOFs"
            << std::setw(14) << "eta" << std::setw(9) << "Rate s" << std::endl;

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
  // One VTK file per refinement level, the tables below are the interesting output
  chemfem::SetVerbose(false);

  GenericEstimator Estimator = GenericEstimator::Residual(f);

  // ------------------------------------------------------------------ uniform
  std::vector<Row> uniform;

  {
    LShapeMesh mesh(4);

    const int nr_levels = 4;

    for(int level=0; level<nr_levels; ++level)
      {
        LagrangeElement element(Degree);
        FESpace Space(mesh, element, WholeBoundary);

        FEFunction Sol = Solve(Space);

        Row row;
        row.cells = mesh.NrCells();
        row.dofs  = Space.NrFreeDof();
        row.eta   = sqrt(Sum(Estimator.Assemble(Sol)));
        uniform.push_back(row);

        if(level+1 < nr_levels)
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

    // Refine until the adaptive mesh has overtaken the uniform reference, so that
    // both tables can be compared at nearly the same number of degrees of freedom
    const size_t target_dofs = uniform.back().dofs;
    const int max_iter = 60;

    for(int level=0; level<max_iter; ++level)
      {
        LagrangeElement element(Degree);
        FESpace Space(mesh, element, WholeBoundary);

        FEFunction Sol = Solve(Space);

        Vector Indicators = Estimator.Assemble(Sol);

        Row row;
        row.cells = mesh.NrCells();
        row.dofs  = Space.NrFreeDof();
        row.eta   = sqrt(Sum(Indicators));
        adaptive.push_back(row);

        Sol.WriteVtk("adaptive_" + std::to_string(level) + ".vtk");

        grading = Grading(mesh);

        if(row.dofs > target_dofs)
          break;

        mesh.Refine(MarkCells(Indicators, 0.6));

        if(!mesh.Check())
          {
            std::cerr << "ERROR: the mesh is broken after refinement step "
                      << level+1 << ".\n";
            return 1;
          }
      }
  }

  PrintTable("Uniform refinement", uniform);
  PrintTable("Adaptive refinement (Doerfler, theta = 0.6)", adaptive);

  // The adaptive row closest in size to the finest uniform one
  const Row& reference = uniform.back();

  size_t best = 0;
  for(size_t i=1; i<adaptive.size(); ++i)
    if(std::labs(long(adaptive[i].dofs) - long(reference.dofs))
       < std::labs(long(adaptive[best].dofs) - long(reference.dofs)))
      best = i;

  // Averaged over the last few rows, a single step is too noisy
  const size_t span = adaptive.size() > 5 ? 5 : adaptive.size()-1;
  const double rate = Rate(adaptive[adaptive.size()-1-span], adaptive.back());
  const double gain = reference.eta / adaptive[best].eta;

  std::cout << "\n" << std::fixed << std::setprecision(3)
            << "Rate adaptive                 : " << rate
            << "   (expected " << OptimalRate << ", optimal for P" << Degree << ")\n"
            << "eta uniform / eta adaptive    : " << gain
            << "   at " << reference.dofs << " against "
            << adaptive[best].dofs << " DOFs\n"
            << "Mesh grading h_max / h_min    : " << grading << std::endl;

  bool ok = true;

  if(adaptive.back().dofs <= reference.dofs)
    {
      std::cerr << "ERROR: the adaptive loop did not reach the size of the uniform "
                << "reference, the two tables are not comparable.\n";
      ok = false;
    }

  // A missing Laplacian only degrades the rate, so the tolerance stays tight
  if(std::fabs(rate - OptimalRate) > 0.1*Degree)
    {
      std::cerr << "ERROR: the adaptive rate is not the expected "
                << OptimalRate << ".\n";
      ok = false;
    }

  if(gain < 2.)
    {
      std::cerr << "ERROR: adaptive refinement does not pay off against uniform.\n";
      ok = false;
    }

  // Uniform refinement keeps h_max/h_min at exactly 1, so anything well above that
  // proves the refinement concentrated
  if(grading < 8.)
    {
      std::cerr << "ERROR: the mesh is not graded, the refinement spread out.\n";
      ok = false;
    }

  if(!ok) return 1;

  std::cout << "\nAdaptivityTest was successful.\n";
  return 0;
}
