#include "fem/FEExpression.h"

namespace chemfem{
  namespace fem{

    FEExpression::FEExpression(ExpressionType type, ScalarFunction Coeff)
      : type(type), Coeff(Coeff) {}
    
  }
}
