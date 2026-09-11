#ifndef _BILINEAR_FORM_H_
#define _BILINEAR_FORM_H_

#include "linalg/SparseMatrix.h"
#include "linalg/SparseMatrixInserter.h"
#include "fem/FESpace.h"
#include "fem/FEExpression.h"

using chemfem::linalg::SparseMatrix;

namespace chemfem{
  namespace fem{

    /**
     * This class represents a finite element bilinear form. In  linear algebra context
     * this corresponds to a matrix (e.g. stiffness or mass matrix).
     */
    class BilinearForm
    {
    public:

      /**
       * Constructor which initializes an empty bilinear form for a given
       * trial and test space.
       */
      BilinearForm(const FESpace&, const FESpace&);

      /// Adds a Laplace term (\nabla u,\nabla v) to the bilinear form
      void AddLaplaceTerm();

      /**
       * Adds the diffusion term div(a\nabla u) to the bilinear form.
       * The coefficient a:R->R passed to the function is the diffusion coefficient.
       * In case of a(x)=1 this corresponds to the Laplace operator. In this case
       * the function AddLaplaceTerm should be used.
       */
      void AddDiffusionTerm(ScalarFunction);

      /**
       * Adds a convection term (b.grad u, v) to the bilinear form. The convection
       * field b is given as a vector valued coefficient function.
       */
      void AddConvectionTerm(VectorFunction);

      /**
       * Adds a reaction term (c u,v) to the bilinear form. The reaction parameter
       * is given as a coefficient function.
       */
      void AddReactionTerm(ScalarFunction);

      /**
       * Adds the term (c A u, B v), where the operators A and B, each one of VALUE, DX
       * and DY, act on the trial and the test function. The Laplace term, for example,
       * is (DX u, DX v) + (DY u, DY v).
       */
      void AddTerm(ScalarFunction, FEOperator, FEOperator);

      /// Adds the term (A u, B v) with the coefficient 1
      void AddTerm(FEOperator, FEOperator);

      /**
       * Assembles the finite element matrix.
       */
      void Assemble();

      /**
       * Assembles the matrix into a larger one, starting at the given row and column.
       * Used for the blocks of a BlockSystem.
       */
      void Assemble(chemfem::linalg::SparseMatrixInserter&, size_t RowOffset, size_t ColOffset);

      /**
       * Returns the matrix which corresponds to the bilinear form. Before calling this function
       * the assemble routine has to be invoked. Otherwise, an empty matrix is returned.
       */
      SparseMatrix& SystemMatrix();

      const FESpace& GetTrialSpace() const;

      const FESpace& GetTestSpace() const;

    private:
      const FESpace& TrialSpace;
      const FESpace& TestSpace;

      SparseMatrix Matrix;
      Vector DirichletRhs;

      std::vector<FEExpression> Terms;
    };

  };
};

#endif
