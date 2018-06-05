#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include "linalg/Vector.h"

using chemfem::linalg::Vector;

namespace chemfem{
  namespace fem{

    /// Finite element type.
    enum FEType {Lagrange};
    
    /**
     * This class represents a single finite element. This is a virtual class and 
     * one should use some child class.
     */
    class Element{

    public:

      /**
       * Initialize by type and degree.
       */
      Element(FEType, int);
      
      /**
       * Returns the number of local degrees of freedom.
       */
      int NrDof() const;

      /**
       * Returns the type of the finite element.
       */
      FEType Type() const;

      /**
       * Returns the degree of the finite element.
       */
      int Degree() const;
      
      /**
       * Return the function value of the trial functions
       */
      virtual double Value(int, double, double) const = 0;

      /**
       * Returns the gradient of the trial function.
       */ 
      virtual Vector Gradient(int, double, double) const = 0;
      
    protected:
      FEType type;
      int nr_dof;
      int degree;
    };
    
  };
};

#endif
