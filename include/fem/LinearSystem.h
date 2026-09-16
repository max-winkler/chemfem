#ifndef _LINEAR_SYSTEM_H_
#define _LINEAR_SYSTEM_H_

#include "fem/BlockSystem.h"

namespace chemfem{
  namespace fem{

    /**
     * Linear system for a single unknown. It is a BlockSystem with one block, without the
     * indices that carry no information there:
     *
     * \code
     *   LinearSystem S(V);
     *   S.Add(A);                   // a bilinear form, several of them add up
     *   S.Add(F);                   // a linear form on the right hand side
     *   S.SetDirichletValues(g);    // named once, the lifting is taken care of
     *   S.AssembleMatrix();
     *
     *   FEFunction u = S.Extract(S.Solve(S.AssembleRhs()));
     * \endcode
     *
     * The prescribed values are named once here instead of once for the form, once for the
     * lifting of the right hand side and once for the solution. A space without a Dirichlet
     * boundary leaves the system singular, which FixDof or AddMeanValueConstraint repairs.
     *
     * Solve uses UMFPACK. To choose a solver, assemble and work with SystemMatrix() instead.
     */
    class LinearSystem
    {
    public:
      explicit LinearSystem(const FESpace&);

      /// Adds a bilinear form to the system, several forms add up
      void Add(BilinearForm&);

      /// Adds a linear form to the right hand side, several forms add up
      void Add(LinearForm&);

      /**
       * Prescribes the Dirichlet values. They enter the right hand side as the lifting of
       * the matrix columns, and Extract puts them into the solution.
       */
      void SetDirichletValues(const DirichletValues&);

      /**
       * Fixes one degree of freedom to zero, which determines the constant of a problem with
       * natural conditions only. This keeps the matrix sparse, unlike AddMeanValueConstraint,
       * and the solution is normalized afterwards with FEFunction::SubtractMean.
       */
      void FixDof(size_t dof = 0);

      /**
       * Fixes the constant by requiring that the integral of the solution vanishes. The
       * constraint is imposed with a Lagrange multiplier, which couples all DOFs and makes
       * the factorization considerably more expensive than FixDof.
       */
      void AddMeanValueConstraint();

      /// Assembles the matrix from the bilinear forms and the constraints
      SparseMatrix& AssembleMatrix();

      /**
       * Assembles the right hand side from the linear forms and the Dirichlet lifting. It
       * does not depend on the matrix, so with data that does not change in time it can be
       * computed once outside a time loop.
       */
      Vector AssembleRhs();

      /**
       * Adds a vector to a right hand side. Used for contributions that are not an integral
       * over the mesh, e.g. the mass matrix times the solution of the previous time step.
       */
      void AddToRhs(Vector&, const Vector&) const;

      /// The free degrees of freedom, taken from a solution of the system
      Vector FreeDof(const Vector&) const;

      /**
       * Solves the system for the given right hand side with UMFPACK. The factorization is
       * computed on the first call and kept until the matrix is assembled again, so further
       * right hand sides only cost a forward and a backward substitution.
       */
      Vector Solve(const Vector&);

      SparseMatrix& SystemMatrix();

      /// Number of unknowns, without the Lagrange multiplier of a mean value constraint
      size_t NrDof() const;

      /// The FE function of the solution
      FEFunction Extract(const Vector&) const;

    private:
      BlockSystem System;
    };

  };
};

#endif
