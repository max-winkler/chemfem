#include "fem/Element.h"

namespace chemfem{
  namespace fem{

    Element::Element(FEType type, int degree) : type(type), degree(degree) {}
    
    int Element::NrDof()
    {
      return nr_dof;
    }
    
  }
}
