#ifndef _DIRICHLET_VALUES_H_
#define _DIRICHLET_VALUES_H_

#include <vector>

#include "fem/FEExpression.h"
#include "fem/FESpace.h"

namespace chemfem{
  namespace fem{

    /**
     * The values prescribed on the Dirichlet boundary of a space. The space decides which
     * degrees of freedom are fixed, this class says which value each of them takes:
     *
     * \code
     *   FESpace V(mesh, element, IsDirichlet);
     *
     *   DirichletValues g(V);        // zero on the whole Dirichlet boundary
     *   g.Set(inflow, IsInflow);     // a profile on one part of it
     *
     *   A.SetDirichletValues(g);     // the bilinear form computes the lifting A_fd g_d
     *   A.Assemble();
     *
     *   Sol.CreateFunction(A.SystemMatrix().Solve(F.LoadVector() - A.DirichletRhs()), g);
     * \endcode
     *
     * The values are taken in the points of the degrees of freedom, so the element has to
     * have point values as DOFs, as Lagrange and Crouzeix-Raviart do.
     */
    class DirichletValues
    {
    public:
      /// Zero on the whole Dirichlet boundary of the space
      explicit DirichletValues(const FESpace&);

      /**
       * Prescribes the values on the part of the Dirichlet boundary where the indicator is
       * true, on all of it if the indicator is omitted
       */
      void Set(ScalarFunction, BoundaryIndicator = nullptr);

      /// Value prescribed for a degree of freedom, zero for the free ones
      double operator[](size_t) const;

      const FESpace& GetFESpace() const;

    private:
      const FESpace* Space;

      /// One entry per DOF of the space, the free ones stay zero
      std::vector<double> Values;
    };

  };
};

#endif
