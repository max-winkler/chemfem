#ifndef _BILINEAR_FORM_H_
#define _BILINEAR_FORM_H_

#include "linalg/DenseMatrix.h"
#include "linalg/SparseMatrix.h"
#include "linalg/SparseMatrixInserter.h"
#include "fem/FESpace.h"
#include "fem/FEExpression.h"
#include "fem/WeakForm.h"

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
       * Adds a volume integral in the notation of WeakForm.h, e.g.
       *   A.AddVolumeTerm(Dx(u)*Dx(v) + Dy(u)*Dy(v))
       */
      void AddVolumeTerm(const BilinearExpression&);

      /**
       * Adds an integral over the part of the boundary where the indicator is true, over
       * the whole boundary if it is omitted, e.g. the Robin term
       *   A.AddBoundaryTerm(alpha * u*v, RobinPart)
       */
      void AddBoundaryTerm(const BilinearExpression&, BoundaryIndicator = nullptr);

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
      /// Adds the local matrix of a cell to the global one, only the free DOFs are kept
      void InsertLocalMatrix(chemfem::linalg::SparseMatrixInserter&, size_t RowOffset,
                             size_t ColOffset, size_t Cell, const chemfem::linalg::DenseMatrix&);

      struct BoundaryProduct
      {
        BilinearProduct product;
        BoundaryIndicator part;
      };

      const FESpace& TrialSpace;
      const FESpace& TestSpace;

      SparseMatrix Matrix;
      Vector DirichletRhs;

      std::vector<FEExpression> Terms;

      std::vector<BilinearProduct> VolumeProducts;
      std::vector<BoundaryProduct> BoundaryProducts;
    };

  };
};

#endif
