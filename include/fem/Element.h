#ifndef _ELEMENT_H_
#define _ELEMENT_H_

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
      int NrDof();

      /**
       * Return the function value of the trial functions
       */
      virtual double Value(int, double, double) = 0;
      
    protected:
      FEType type;
      int nr_dof;
      int degree;
    };
    
  };
};

#endif
