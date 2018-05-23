#include <iostream>

#include "linalg/SparseMatrix.h"

using namespace chemfem::linalg;

int main()
{
  SparseMatrix identity(IdentityMatrix(5));

  std::cout << " identity = \n" << identity << std::endl;
  
  return 0;
}
