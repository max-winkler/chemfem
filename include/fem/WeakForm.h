#ifndef _WEAK_FORM_H_
#define _WEAK_FORM_H_

#include <utility>
#include <vector>

#include "fem/FEExpression.h"
#include "linalg/Coordinate.h"
#include "linalg/Vector2D.h"

namespace chemfem::mesh{
  class Mesh;
}

namespace chemfem{
  namespace fem{

    class FEFunction;

    /*
     * Notation for the terms of bilinear and linear forms, close to the weak formulation:
     *
     *   TrialFunction u;  TestFunction v;
     *
     *   A.AddVolumeTerm(nu * (Dx(u)*Dx(v) + Dy(u)*Dy(v)) + 1./tau * u*v);
     *   A.AddBoundaryTerm(alpha * u*v, RobinPart);
     *   F.AddVolumeTerm(f * v + 1./tau * Uold * v);
     *
     * u and v are only symbols, Dx and Dy take the partial derivatives. Coefficients are
     * written in front and may be numbers, functions of the coordinates, functors or FE
     * functions on the same mesh. An FE function is only referenced, it has to live as
     * long as the form is assembled.
     */

    /// Operators applied to the trial and test functions
    enum FEOperator {VALUE, DX, DY};

    /// Applies the operator to a basis function, given by its value and its gradient
    inline double ApplyOperator(FEOperator op, double value,
                                const chemfem::linalg::Vector2D& grad)
    {
      switch(op)
        {
        case DX: return grad.x;
        case DY: return grad.y;
        default: return value;
        }
    }

    /// Symbol for the trial function of a bilinear form
    struct TrialFunction {};

    /// Symbol for the test function of a bilinear or linear form
    struct TestFunction {};

    /**
     * Coefficient of a term, the product of a number, a function of the coordinates and
     * an FE function. Each factor is optional.
     */
    class Coefficient
    {
    public:
      Coefficient(double);

      Coefficient(const FEFunction&);

      /// Anything that can be called with a Coordinate and returns a number
      template<typename F, typename = decltype(double(std::declval<const F&>()(
                             std::declval<const chemfem::linalg::Coordinate&>())))>
      Coefficient(const F& f) : factor(1.), function(f), fe(nullptr) {}

      /// Value in the point x, which is (xi,eta) on the reference element of the cell
      double Value(const chemfem::linalg::Coordinate& x, size_t cell,
                   double xi, double eta) const;

      /// False, with an error message, if an FE function in it lives on another mesh
      bool LivesOn(const chemfem::mesh::Mesh&) const;

      friend Coefficient operator*(const Coefficient&, const Coefficient&);

    private:
      double factor;
      ScalarFunction function;
      const FEFunction* fe;
    };

    Coefficient operator*(const Coefficient&, const Coefficient&);

    /// Trial function with an operator and a coefficient, e.g. alpha * Dx(u)
    struct TrialTerm
    {
      TrialTerm(TrialFunction);
      TrialTerm(const Coefficient&, FEOperator);

      Coefficient coeff;
      FEOperator op;
    };

    /// Test function with an operator and a coefficient, e.g. f * v
    struct TestTerm
    {
      TestTerm(TestFunction);
      TestTerm(const Coefficient&, FEOperator);

      Coefficient coeff;
      FEOperator op;
    };

    /// Product of a trial and a test function, each with its operator, and a coefficient
    struct BilinearProduct
    {
      Coefficient coeff;
      FEOperator trial, test;
    };

    /// Sum of such products, the terms of a bilinear form
    struct BilinearExpression
    {
      std::vector<BilinearProduct> products;
    };

    /// Sum of test functions with operators and coefficients, the terms of a linear form
    struct LinearExpression
    {
      LinearExpression(const TestTerm&);

      std::vector<TestTerm> terms;
    };

    TrialTerm Dx(TrialFunction);
    TrialTerm Dy(TrialFunction);
    TestTerm Dx(TestFunction);
    TestTerm Dy(TestFunction);

    TrialTerm operator*(const Coefficient&, const TrialTerm&);
    TestTerm operator*(const Coefficient&, const TestTerm&);

    BilinearExpression operator*(const TrialTerm&, const TestTerm&);
    BilinearExpression operator*(const TestTerm&, const TrialTerm&);

    BilinearExpression operator*(const Coefficient&, const BilinearExpression&);
    LinearExpression operator*(const Coefficient&, const LinearExpression&);

    BilinearExpression operator+(const BilinearExpression&, const BilinearExpression&);
    LinearExpression operator+(const LinearExpression&, const LinearExpression&);

  };
};

#endif
