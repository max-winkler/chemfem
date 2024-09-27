#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "fem/LinearForm.h"
#include "fem/BilinearForm.h"
#include "fem/LagrangeElement.h"
#include "fem/FEFunction.h"
#include "fem/ErrorEstimator.h"
#include "mesh/LShapeMesh.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

double f(double x, double y)
{
  return exp(-pow((x+0.5) / 0.2, 2.) - pow((y-0.5)/0.2, 2.));
}


int main()
{
  LShapeMesh mesh(10);

  const int max_iter = 8;

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
      std::cout << "Error of equation system: " << Res.Norm() << std::endl;
      
      FEFunction Sol(Space);
      Sol.CreateFunction(X);
      Sol.WriteVtk("solution.vtk");

      // ESTIMATE
      ErrorEstimator Estimator(Sol);
      Estimator.AddVolumeResidualTerm(f);

      Vector Errors = Estimator.Assemble();

      // MARK
      double Threshold = 0.5;
      
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
    }
}
