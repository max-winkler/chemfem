#include "fem/FESpace.h"

namespace chemfem{
  namespace fem{

    FESpace::FESpace(const Mesh& mesh, FEType type, int degree)
      : mesh(mesh), type(type), degree(degree) {}
    
  }
}
