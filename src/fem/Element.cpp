#include "fem/Element.h"

namespace chemfem{
  namespace fem{

    Element::Element(FEType type, int degree) : type(type), degree(degree) {}
    
    int Element::NrDof() const
    {
      return nr_dof;
    }

    FEType Element::Type() const
    {
      return type;
    }
    
  }
}
