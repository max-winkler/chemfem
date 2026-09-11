#include "fem/FEExpression.h"

namespace chemfem{
  namespace fem{

    FEExpression::FEExpression(ExpressionType Type, ScalarFunction Coeff)
      : Type(Type), Coeff(Coeff), TrialOp(VALUE), TestOp(VALUE) {}

    FEExpression::FEExpression(ExpressionType Type, VectorFunction VecCoeff)
      : Type(Type), VecCoeff(VecCoeff), TrialOp(VALUE), TestOp(VALUE) {}

    FEExpression::FEExpression(ExpressionType Type)
      : Type(Type), Coeff(nullptr), TrialOp(VALUE), TestOp(VALUE) {}

    FEExpression::FEExpression(ScalarFunction Coeff, FEOperator TrialOp, FEOperator TestOp)
      : Type(GENERAL), Coeff(Coeff), TrialOp(TrialOp), TestOp(TestOp) {}

    ExpressionType FEExpression::GetType() const
    {
      return Type;
    }

    double FEExpression::EvalCoeff(const chemfem::linalg::Coordinate& p) const
    {
      return Coeff(p);
    }

    chemfem::linalg::Vector2D FEExpression::EvalVectorCoeff(const chemfem::linalg::Coordinate& p) const
    {
      return VecCoeff(p);
    }
  }
}
