#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "fem/BlockSystem.h"
#include "fem/ProductElement.h"
#include "fem/LagrangeElement.h"
#include "mesh/UnitSquareMesh.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// A system may be assembled thousands of times, once per step of a Newton iteration. Every
// assembly replaces the matrix, so SparseMatrix has to free what it held before. It did not,
// which nobody noticed while every test assembled once, and which killed a Navier-Stokes run
// after roughly 1900 assemblies.

/// (grad u, grad v) for vector valued functions
double Stiffness(const VectorValues& u, const VectorValues& v)
{
  return ddot(u.gradient, v.gradient);
}

/// -(div u, q)
double Divergence(const VectorValues& u, const PointValues& q)
{
  return -u.divergence * q.value;
}

/// Resident memory of this process in MB, negative if it cannot be read
double ResidentMB()
{
  std::ifstream status("/proc/self/status");
  std::string word;

  while(status >> word)
    if(word == "VmRSS:")
      {
        double kb = 0.;
        status >> kb;
        return kb/1024.;
      }

  return -1.;
}

int main()
{
  Mesh mesh = UnitSquareMesh(3);
  mesh.RefineUniform();
  mesh.RefineUniform();

  LagrangeElement P2(2), P1(1);
  ProductElement Velocity(P2, 2);

  FESpace V(mesh, Velocity, WholeBoundary), Q(mesh, P1);

  BilinearForm A(V, V), B(V, Q);
  A.AddVolumeTerm(Stiffness);
  B.AddVolumeTerm(Divergence);

  BlockSystem S({V, Q});
  S.AddBlock(0, 0, A);
  S.AddBlock(1, 0, B);
  S.AddTransposedBlock(0, 1, B);

  const int repeats = 60;

  // The first assemblies still grow the allocator, the slope is taken after them
  const int warmup = 10;

  for(int r=0; r<warmup; ++r)
    S.AssembleMatrix();

  const double before = ResidentMB();

  if(before < 0.)
    {
      std::cout << "LeakTest needs /proc/self/status, which this system does not offer. "
                << "Nothing was checked.\n";
      return 0;
    }

  for(int r=warmup; r<repeats; ++r)
    S.AssembleMatrix();

  const double after = ResidentMB();
  const double per_assembly = (after - before)/(repeats - warmup);

  std::cout << mesh.NrCells() << " cells, " << S.NrDof() << " unknowns, "
            << repeats << " assemblies\n"
            << std::fixed << std::setprecision(1)
            << "Resident memory " << before << " MB before, " << after << " MB after, "
            << std::setprecision(3) << per_assembly << " MB per assembly" << std::endl;

  // Without the destructor and the assignment operator of SparseMatrix every assembly used
  // to leak its three arrays, which is far above this bound on this mesh
  if(per_assembly > 0.05)
    {
      std::cerr << "ERROR: repeated assembly leaks " << per_assembly
                << " MB each time. SparseMatrix has to free its arrays.\n";
      return 1;
    }

  std::cout << "\nLeakTest was successful.\n";
  return 0;
}
