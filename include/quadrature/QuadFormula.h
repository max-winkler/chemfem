#ifndef _QUAD_FORMULA_H_
#define _QUAD_FORMULA_H_

#include "linalg/Vector.h"

using chemfem::linalg::Vector;

namespace chemfem{
  namespace quadrature{

    enum QUAD_FORMULA {MIDPOINT, VERTEX};
    
    /**
     * This class stores all information on a certain quadratue formula.      
     */
    class QuadratureFormula
    {
    public:
      /**
       * Initialize a certain quadrature formula. All possible choices can be found in the 
       * enumeration QUAD_FORMULA.
       */
      QuadratureFormula(QUAD_FORMULA);

      /**
       * Find a suitable quadrature formula automatically which is exact for polynomials 
       * of a certain degree passed as input argument to this function.
       */
      QuadratureFormula(int);
      
      /**
       * Returns the number of quadrature points.
       */ 
      size_t NrQuadPoints();

      /**
       * Sets the pointers to the quadrature formula data.
       * The first argument are the weights, the second and third argument are the the 
       * coordinates of quadrature points. For a 1D quadrature formula the third argument is ignored.
       */
      void FormulaData(Vector&, Vector&, Vector&);
      
    private:
      int Points;
      Vector Weights;
      Vector Xi, Eta;
    };
  };
};

#endif
