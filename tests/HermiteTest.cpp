#include "fem/BilinearForm.h"
#include "fem/FEExpression.h"
#include "fem/FEFunction.h"
#include "fem/HermiteElement.h"
#include "fem/LinearForm.h"
#include "fem/LinearSystem.h"
#include "fem/VtkOutput.h"
#include "mesh/Cell.h"
#include "mesh/UnitSquareMesh.h"
#include "fem/FESpace.h"

using namespace chemfem::linalg;
using namespace chemfem::fem;
using namespace chemfem::mesh;

int main() {
  chemfem::mesh::UnitSquareMesh mesh(10);

  chemfem::fem::HermiteElement element;

  FESpace V(mesh, element, WholeBoundary);

  BilinearForm a(V, V);
  a.AddLaplaceTerm();

  LinearForm F(V);
  F.AddVolumeForce([](const Coordinate& v) { return 1.; });

  // Homogeneous conditions need no DirichletValues: the constrained DOFs stay zero
  LinearSystem S(V);
  S.AddLhs(a);
  S.AddRhs(F);
  S.AssembleMatrix();

  FEFunction y = S.Extract(S.Solve(S.AssembleRhs()));

  VtkOutput output(mesh);
  output.AddScalar("y", y);
  output.Write("poisson_hermite.vtk");
  
  return 0;  
}  
