#include <iostream>

#include "linalg/SparseMatrix.h"
#include "linalg/SparseMatrixInserter.h"

using namespace chemfem::linalg;

int main()
{
  SparseMatrix identity(IdentityMatrix(5));

  std::cout << " identity = \n" << identity << std::endl;

  SparseMatrix custom(3, 5);
  SparseMatrixInserter ins(custom);

  ins.Insert(0,0,1);
  ins.Insert(0,1,2);
  ins.Insert(1,2,3);
  ins.Insert(2,2,4);
  ins.Insert(2,4,5);
  ins.Build();

  std::cout << " curstom = \n" << custom << std::endl;
  
  return 0;
}
