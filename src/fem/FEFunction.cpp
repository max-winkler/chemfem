#include "fem/FEFunction.h"

namespace chemfem{
  namespace fem{
    
    FEFunction::FEFunction(const FEFunction& other)
      : Space(other.Space), Data(other.Data) {}
    
  }
}
