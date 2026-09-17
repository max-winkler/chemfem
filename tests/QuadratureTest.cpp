#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

#include "quadrature/QuadFormula.h"

using namespace chemfem::quadrature;
using chemfem::linalg::Vector;

// Every formula on the triangle has to integrate the monomials up to its degree exactly.
// The integral of xi^p eta^q over the reference triangle is p! q! / (p+q+2)!, which gives a
// sharp check without any finite element around it.

double Factorial(int n)
{
  double f = 1.;
  for(int k=2; k<=n; ++k)
    f *= k;
  return f;
}

double Exact(int p, int q)
{
  return Factorial(p)*Factorial(q)/Factorial(p + q + 2);
}

/// Largest relative error of the formula on all monomials of the given total degree
double WorstError(QUAD_FORMULA formula, int total_degree)
{
  QuadratureFormula Quad(formula);

  Vector Weights, Xi, Eta;
  Quad.FormulaData(Weights, Xi, Eta);

  double worst = 0.;

  for(int p=0; p<=total_degree; ++p)
    {
      const int q = total_degree - p;

      double sum = 0.;
      for(size_t k=0; k<Weights.size(); ++k)
        sum += Weights[k] * std::pow(Xi[k], p) * std::pow(Eta[k], q);

      worst = std::max(worst, std::fabs(sum - Exact(p, q))/Exact(p, q));
    }

  return worst;
}

bool Check(const std::string& name, QUAD_FORMULA formula, int degree, double tolerance)
{
  QuadratureFormula Quad(formula);

  std::cout << std::setw(14) << name << std::setw(10) << Quad.NrQuadPoints()
            << std::setw(10) << degree;

  bool ok = true;

  // Exact up to its degree
  double worst = 0.;
  for(int d=0; d<=degree; ++d)
    worst = std::max(worst, WorstError(formula, d));

  std::cout << std::scientific << std::setprecision(2) << std::setw(14) << worst;

  if(worst > tolerance)
    {
      std::cerr << "\nERROR: " << name << " is not exact up to degree " << degree
                << " within " << tolerance << ".\n";
      ok = false;
    }

  // One degree higher it must not be exact any more, otherwise the stated degree is too low
  // and the formula is more expensive than it needs to be
  const double beyond = WorstError(formula, degree+1);

  std::cout << std::setw(14) << beyond << std::endl;

  if(beyond < 1.e-13)
    {
      std::cerr << "\nWarning: " << name << " is exact beyond degree " << degree
                << " as well.\n";
    }

  return ok;
}

int main()
{
  std::cout << std::setw(14) << "formula" << std::setw(10) << "points"
            << std::setw(10) << "degree" << std::setw(14) << "worst error"
            << std::setw(14) << "next degree" << "\n"
            << std::string(62, '=') << std::endl;

  bool ok = true;

  ok &= Check("MIDPOINT", MIDPOINT, 1, 1.e-13);
  ok &= Check("GAUSS_EDGE", GAUSS_EDGE, 2, 1.e-13);
  ok &= Check("GAUSS_5", GAUSS_5, 5, 1.e-13);

  // The Xi values of GAUSS_7 are tabulated with ten digits only, see QuadFormula.cpp, while
  // its Eta values carry sixteen. That truncation, not the formula itself, limits it to
  // about 1e-9. Recomputing the sixteen constants would lift it to round-off.
  ok &= Check("GAUSS_7", GAUSS_7, 7, 1.e-8);

  // The degree constructor has to pick the cheapest formula that is exact
  const size_t expected[8] = {1, 1, 3, 7, 7, 7, 16, 16};

  std::cout << "\nFormula chosen by degree:";

  for(int degree=0; degree<=7; ++degree)
    {
      QuadratureFormula Quad(degree);
      const size_t points = Quad.NrQuadPoints();

      std::cout << " " << degree << "->" << points;

      if(points != expected[degree])
        {
          std::cerr << "\nERROR: degree " << degree << " gave a formula with " << points
                    << " points, expected " << expected[degree] << ".\n";
          ok = false;
        }
    }

  std::cout << std::endl;

  if(!ok)
    return 1;

  std::cout << "\nQuadratureTest was successful.\n";
  return 0;
}
