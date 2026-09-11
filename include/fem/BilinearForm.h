#ifndef _BILINEAR_FORM_H_
#define _BILINEAR_FORM_H_

#include <functional>

#include "linalg/DenseMatrix.h"
#include "linalg/SparseMatrix.h"
#include "linalg/SparseMatrixInserter.h"
#include "fem/FESpace.h"
#include "fem/FEExpression.h"
#include "fem/PointValues.h"

using chemfem::linalg::SparseMatrix;

namespace chemfem{
  namespace fem{

    /**
     * This class represents a finite element bilinear form. In  linear algebra context
     * this corresponds to a matrix (e.g. stiffness or mass matrix).
     *
     * Terms beyond the predefined ones are given as functors, the same way as for the
     * GenericEstimator. The functor gets the values of the trial function u and the
     * test function v in a quadrature point:
     *
     * \code
     *   struct Mass
     *   {
     *     double operator()(const QuadPoint&, const CellGeometry&,
     *                       const PointValues& u, const PointValues& v) const
     *     {
     *       return u.value * v.value;
     *     }
     *   };
     *
     *   A.AddVolumeTerm(Mass());
     * \endcode
     */
    class BilinearForm
    {
    public:

      /// Integrand of a volume term for the trial function u and the test function v
      typedef std::function<double(const QuadPoint& p, const CellGeometry& cell,
                                   const PointValues& u, const PointValues& v)> VolumeIntegrand;

      /// Integrand of a boundary term, evaluated on a boundary edge of the cell
      typedef std::function<double(const QuadPoint& p, const CellGeometry& cell,
                                   const EdgeGeometry& edge,
                                   const PointValues& u, const PointValues& v)> BoundaryIntegrand;

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

      /// Adds the integral of the functor over all cells
      void AddVolumeTerm(VolumeIntegrand);

      /**
       * Adds the integral of the functor over the part of the boundary where the
       * indicator is true, over the whole boundary if it is omitted
       */
      void AddBoundaryTerm(BoundaryIntegrand, BoundaryIndicator = nullptr);

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

      struct BoundaryTerm
      {
        BoundaryIntegrand integrand;
        BoundaryIndicator part;
      };

      const FESpace& TrialSpace;
      const FESpace& TestSpace;

      SparseMatrix Matrix;
      Vector DirichletRhs;

      std::vector<FEExpression> Terms;

      std::vector<VolumeIntegrand> VolumeTerms;
      std::vector<BoundaryTerm> BoundaryTerms;
    };

  };
};

#endif
