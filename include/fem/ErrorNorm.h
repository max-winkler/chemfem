#ifndef _ERROR_NORM_H_
#define _ERROR_NORM_H_

#include <iostream>

#include "fem/FEFunction.h"
#include "fem/FEExpression.h"

namespace chemfem{
  namespace fem{

    enum Norm{L2, H1, H1_SEMI};
    
    /**
     * This class provides the functionality to compute the error between 
     * the finite element solution and the exact solution of a problem.
     */
    class ErrorNorm
    {
    public:

      /**
       * Provide an exact solution of the problem.
       */
      void SetExactValue(ScalarFunction);

      /**
       * Provide the first derivatives of the exact solution.
       */
      void SetExactGradient(VectorFunction);

      /**
       * Set the finite element solution which should be compared to the exact solution.
       */
      void SetFEFunction(FEFunction&);
      
      /**
       * Compute the error Norm. The argument specifies the norm in which the error
       * should be evaluated.
       */
      double Compute(Norm) const;
      
    private:
      
      ScalarFunction Value;
      VectorFunction Gradient;
      const FEFunction* FESolution;
    };

  };
};

#endif
