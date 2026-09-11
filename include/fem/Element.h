#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include "linalg/Matrix2D.h"
#include "linalg/Vector2D.h"

using chemfem::linalg::Matrix2D;
using chemfem::linalg::Vector2D;

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
      virtual Vector2D Gradient(int, double, double) const = 0;

      /**
       * Returns the Hessian of the trial function on the reference element. It is
       * symmetric, so only one of the two off-diagonal entries carries information.
       */
      virtual Matrix2D Hessian(int, double, double) const = 0;
      
    protected:
      FEType type;
      int nr_dof;
      int degree;
    };
    
  };
};

#endif
