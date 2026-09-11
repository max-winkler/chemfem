#include "fem/WeakForm.h"

#include <iostream>

#include "fem/FEFunction.h"

using chemfem::linalg::Coordinate;

namespace chemfem{
  namespace fem{

    Coefficient::Coefficient(double value) : factor(value), fe(nullptr) {}

    Coefficient::Coefficient(const FEFunction& u) : factor(1.), fe(&u) {}

    double Coefficient::Value(const Coordinate& x, size_t cell, double xi, double eta) const
    {
      double value = factor;

      if(function)
        value *= function(x);

      if(fe)
        {
          const FESpace& Space = fe->GetFESpace();

          double u = 0.;
          for(size_t k=0; k<Space.NrLocalDof(); ++k)
            u += (*fe)[Space.GetGlobalIndex(cell, k)] * Space.RefElement().Value(k, xi, eta);

          value *= u;
        }

      return value;
    }

    bool Coefficient::LivesOn(const chemfem::mesh::Mesh& mesh) const
    {
      if(!fe || &fe->GetFESpace().GetMesh() == &mesh)
        return true;

      std::cerr << "Error: The FE function in a coefficient has to live on the mesh of the "
                << "form.\n";
      return false;
    }

    Coefficient operator*(const Coefficient& a, const Coefficient& b)
    {
      Coefficient c(a.factor * b.factor);

      if(a.function && b.function)
        {
          const ScalarFunction fa = a.function, fb = b.function;
          c.function = [fa, fb](const Coordinate& x) { return fa(x)*fb(x); };
        }
      else
        c.function = a.function ? a.function : b.function;

      if(a.fe && b.fe)
        std::cerr << "Error: A coefficient can contain only one FE function.\n";

      c.fe = a.fe ? a.fe : b.fe;

      return c;
    }

    TrialTerm::TrialTerm(TrialFunction) : coeff(1.), op(VALUE) {}

    TrialTerm::TrialTerm(const Coefficient& coeff, FEOperator op) : coeff(coeff), op(op) {}

    TestTerm::TestTerm(TestFunction) : coeff(1.), op(VALUE) {}

    TestTerm::TestTerm(const Coefficient& coeff, FEOperator op) : coeff(coeff), op(op) {}

    LinearExpression::LinearExpression(const TestTerm& term) : terms(1, term) {}

    TrialTerm Dx(TrialFunction) { return TrialTerm(1., DX); }
    TrialTerm Dy(TrialFunction) { return TrialTerm(1., DY); }
    TestTerm Dx(TestFunction) { return TestTerm(1., DX); }
    TestTerm Dy(TestFunction) { return TestTerm(1., DY); }

    TrialTerm operator*(const Coefficient& c, const TrialTerm& u)
    {
      return TrialTerm(c * u.coeff, u.op);
    }

    TestTerm operator*(const Coefficient& c, const TestTerm& v)
    {
      return TestTerm(c * v.coeff, v.op);
    }

    BilinearExpression operator*(const TrialTerm& u, const TestTerm& v)
    {
      BilinearExpression e;
      e.products.push_back(BilinearProduct{u.coeff * v.coeff, u.op, v.op});
      return e;
    }

    BilinearExpression operator*(const TestTerm& v, const TrialTerm& u)
    {
      return u * v;
    }

    BilinearExpression operator*(const Coefficient& c, const BilinearExpression& e)
    {
      BilinearExpression scaled = e;
      for(size_t p=0; p<scaled.products.size(); ++p)
        scaled.products[p].coeff = c * scaled.products[p].coeff;
      return scaled;
    }

    LinearExpression operator*(const Coefficient& c, const LinearExpression& e)
    {
      LinearExpression scaled = e;
      for(size_t t=0; t<scaled.terms.size(); ++t)
        scaled.terms[t].coeff = c * scaled.terms[t].coeff;
      return scaled;
    }

    BilinearExpression operator+(const BilinearExpression& a, const BilinearExpression& b)
    {
      BilinearExpression sum = a;
      sum.products.insert(sum.products.end(), b.products.begin(), b.products.end());
      return sum;
    }

    LinearExpression operator+(const LinearExpression& a, const LinearExpression& b)
    {
      LinearExpression sum = a;
      sum.terms.insert(sum.terms.end(), b.terms.begin(), b.terms.end());
      return sum;
    }

  }
}
