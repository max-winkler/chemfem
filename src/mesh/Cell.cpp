#include "mesh/Cell.h"

namespace chemfem{
  namespace mesh{

    Cell::Cell(size_t n1, size_t n2, size_t n3)
    {
      LocIndex[0] = n1;
      LocIndex[1] = n2;
      LocIndex[2] = n3;
    }
    
  };
};
