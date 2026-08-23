#ifndef _FE_EXPRESSION_H_
#define _FE_EXPRESSION_H_

#include <functional>
#include <vector>

#include "linalg/Coordinate.h"
#include "linalg/Vector2D.h"

namespace chemfem{
  namespace fem{

    /**
     * Coefficients, right hand sides and exact solutions are passed as std::function,
     * so that they may carry state: a plain function, a lambda with captures, a
     * functor or a bound member function are all accepted. Plain function pointers
     * convert implicitly, so &f keeps working wherever it worked before.
     */
    typedef std::function<double(const chemfem::linalg::Coordinate&)> ScalarFunction;
    typedef std::function<chemfem::linalg::Vector2D(const chemfem::linalg::Coordinate&)> VectorFunction;
    
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
      FEExpression(ExpressionType); 
      FEExpression(ExpressionType, ScalarFunction);
      FEExpression(ExpressionType, VectorFunction);

      // \todo These functions are added to avoid these friend declarations. Use these functions in BilinearForm and LinearForm too.
      ExpressionType GetType() const;
      double EvalCoeff(const chemfem::linalg::Coordinate&) const;
      chemfem::linalg::Vector2D EvalVectorCoeff(const chemfem::linalg::Coordinate&) const;
        
    private:
      ExpressionType Type;
      ScalarFunction Coeff;
      /// Used instead of Coeff by the terms whose coefficient is vector valued
      VectorFunction VecCoeff;
    };
    
  };
};

#endif
