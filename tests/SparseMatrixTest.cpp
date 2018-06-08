#include <iostream>

#include "linalg/SparseMatrix.h"
#include "linalg/SparseMatrixInserter.h"

using namespace chemfem::linalg;

int main()
{
  SparseMatrix identity(IdentityMatrix(5));

  std::cout << " identity = \n" << identity << std::endl;

  SparseMatrix A(3, 5);
  SparseMatrixInserter ins(A);

  ins.Insert(0,0,1);
  ins.Insert(0,1,2);
  ins.Insert(1,2,3);
  ins.Insert(2,2,4);
  ins.Insert(2,4,5);
  ins.Build();

  std::cout << " A = \n" << A << std::endl;

  Vector x(5);
  x[0] = 1.;
  x[2] = 2.;
  x[4] = 3.;

  std::cout << "x = \n" << x << std::endl;
  
  Vector y(A*x);

  std::cout << "A*x = \n" << y << std::endl;
  
  return 0;
}
