#ifndef _LAGRANGE_ELEMENT_H_
#define _LAGRANGE_ELEMENT_H_

#include "fem/Element.h"

namespace chemfem{
  namespace fem{

    /**
     * This class represents a Lagrange element of arbitrary order.
     */
    class LagrangeElement : public Element
    {
    public: 
      /**
       * Constructor which initializes the finite element by its degree.
       */
      LagrangeElement(int);

      /**
       * Returns the function value of the ansatz function on the reference element.
       */
      double Value(int, double, double) const;

      Vector Gradient(int, double, double) const;

    };    
  };
};


#endif
