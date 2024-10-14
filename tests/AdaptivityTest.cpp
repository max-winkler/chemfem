#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "fem/LinearForm.h"
#include "fem/BilinearForm.h"
#include "fem/LagrangeElement.h"
#include "fem/FEFunction.h"
#include "fem/ErrorEstimator.h"
#include "mesh/UnitSquareMesh.h"
#include "mesh/LShapeMesh.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

double f(double x, double y)
{
  return x+y;
  return exp(-pow((x+0.5) / 0.2, 2.) - pow((y-0.5)/0.2, 2.));
}


int main()
{
  LShapeMesh mesh(3);
  std::cout << mesh << std::endl;
  const int max_iter = 20;

  for(int iter=0; iter<max_iter; ++iter)
    {
      std::cout << " Computation on level " << iter+1 << std::endl;
      std::cout << "======================================\n";
      
      std::cout << "Nr of nodes : " << mesh.NrNodes() << std::endl;
      std::cout << "Nr of cells : " << mesh.NrCells() << std::endl;

      // SOLVE
      LagrangeElement element(1);
      FESpace Space(mesh, element);

      BilinearForm Laplace(Space, Space);
      Laplace.AddLaplaceTerm();
      Laplace.Assemble();

      SparseMatrix& Matrix = Laplace.SystemMatrix();
      
      LinearForm F(Space);
      F.AddVolumeForce(f);
      F.Assemble();
  
      Vector& Vec = F.LoadVector();

      Vector X(Matrix.Solve(Vec));
      
      Vector Res(Matrix*X - Vec);
      
      FEFunction Sol(Space);
      Sol.CreateFunction(X);
      Sol.WriteVtk("solution_" + std::to_string(iter) + ".vtk");

      // ESTIMATE
      ErrorEstimator Estimator(Sol);
      Estimator.AddVolumeResidualTerm(f);
      Estimator.AddEdgeJumpTerm();
        
      Vector Errors = Estimator.Assemble();

      // MARK
      double Threshold = 0.8;
      
      std::vector<bool> Marker(mesh.NrCells(), false);

      double MaxError = 0.;
      for(auto Error = Errors.begin(); Error != Errors.end(); ++Error)
        if(*Error > MaxError) MaxError = *Error;

      size_t CellInd = 0;
      for(auto Error = Errors.begin(); Error != Errors.end(); ++Error, ++CellInd)
        if(*Error > Threshold*MaxError)
          Marker[CellInd] = true;

      // REFINE
      mesh.Refine(Marker);
      if(!mesh.Check())
        std::cerr << "ERROR: There is a problem with the mesh datastructure.\n";
    }
}
