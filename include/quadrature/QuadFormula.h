#ifndef _QUAD_FORMULA_H_
#define _QUAD_FORMULA_H_

#include "linalg/Vector.h"

using chemfem::linalg::Vector;

namespace chemfem{
  namespace quadrature{

    /**
     * The formulas on the triangle and the degree of polynomials they integrate exactly:
     * MIDPOINT 1 point, degree 1; VERTEX 3 points, degree 1; GAUSS_EDGE 3 points, degree 2;
     * GAUSS_5 7 points, degree 5; GAUSS_7 16 points, degree 7. LINE_GAUSS_5 has 5 points on
     * [0,1] and is exact up to degree 9.
     */
    enum QUAD_FORMULA {MIDPOINT, VERTEX, GAUSS_EDGE, GAUSS_5, GAUSS_7, LINE_GAUSS_5};

    /**
     * This class stores all information on a certain quadratue formula. The formulas on
     * the reference triangle have the weight sum 1/2, its area, and the line formulas
     * live on [0,1] with the weight sum 1.
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
