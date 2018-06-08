#ifndef _FE_EXPRESSION_H_
#define _FE_EXPRESSION_H_

#include <vector>

namespace chemfem{
  namespace fem{

    typedef double (*ScalarFunction)(double, double);
    
    enum ExpressionType {SECOND_ORDER, FIRST_ORDER, ZERO_ORDER, VOLUME_FORCE, NEUMANN_BC};

    /**
     * This class is used to store a single term in a partial differential equation.
     * Here, we distinguish among 2nd, 1st and zero-order terms. 
     * Coefficients belonging to the terms are stored in the class as well.
     */
    class FEExpression
    {
      friend class BilinearForm;
      friend class LinearForm;
      
    public:
      FEExpression(ExpressionType, ScalarFunction); 
    private:
      ExpressionType Type;
      ScalarFunction Coeff;
    };
    
  };
};

#endif
