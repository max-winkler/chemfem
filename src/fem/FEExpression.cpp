#include "fem/FEExpression.h"

namespace chemfem{
  namespace fem{

    FEExpression::FEExpression(ExpressionType Type, ScalarFunction Coeff)
      : Type(Type), Coeff(Coeff) {}

    FEExpression::FEExpression(ExpressionType Type)
      : Type(Type), Coeff(nullptr) {}
    
  }
}
