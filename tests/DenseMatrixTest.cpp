#include <iostream>

#include "linalg/DenseMatrix.h"

using namespace chemfem::linalg;

int main()
{
  DenseMatrix A(2,2);
  const double data[] = {1.0, 2.0, 3.0, 4.0};
  A.Init(data);

  Vector b(2);
  b[0] = 1.;
  b[1] = -1.;

  Vector c(2);
  c = A*b;

  std::cout << "A*b is\n";
  
  for(Vector::const_iterator it = c.begin(); it != c.end(); ++it)
    std::cout << " " << *it << std::endl;

  Vector d(2);
  d[0] = 0;
  d[1] = 1;

  Vector e(2);
  e = c+d;

  double norm = e.Norm();
  if(std::abs(norm) > 1.e-8)
    {
      std::cerr << "The computation for the dense matrices was wrong.\n";
      return -1;
    }
  else
    std::cout << "DenseMatrixTest was successful.\n";
  
  return 0;
}
