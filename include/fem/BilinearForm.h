#ifndef _BILINEAR_FORM_H_
#define _BILINEAR_FORM_H_

#include <functional>

#include "linalg/DenseMatrix.h"
#include "linalg/SparseMatrix.h"
#include "linalg/SparseMatrixInserter.h"
#include "fem/FESpace.h"
#include "fem/FEExpression.h"
#include "fem/PointValues.h"
#include "fem/DirichletValues.h"

using chemfem::linalg::SparseMatrix;

namespace chemfem{
  namespace fem{

    /**
     * This class represents a finite element bilinear form. In  linear algebra context
     * this corresponds to a matrix (e.g. stiffness or mass matrix).
     *
     * Terms beyond the predefined ones are given as an integrand, which gets the values of
     * the trial function u and the test function v in a quadrature point. A plain
     * function, a functor or a lambda may serve as integrand:
     *
     * \code
     *   double Mass(const PointValues& u, const PointValues& v) { return u.value * v.value; }
     *
     *   A.AddVolumeTerm(Mass);
     * \endcode
     *
     * An integrand that depends on the position, e.g. through a coefficient or an FE
     * function, gets the quadrature point as first argument.
     */
    class BilinearForm
    {
    public:

      /// Integrand for the trial function u and the test function v
      typedef std::function<double(const PointValues& u, const PointValues& v)> Integrand;

      /// Integrand of a volume term that depends on the quadrature point as well
      typedef std::function<double(const QuadPoint& p,
                                   const PointValues& u, const PointValues& v)> PointIntegrand;

      /// Integrand of a boundary term that depends on the quadrature point or the edge
      typedef std::function<double(const QuadPoint& p, const EdgeGeometry& edge,
                                   const PointValues& u, const PointValues& v)> BoundaryIntegrand;

      /// Integrand on a vector valued trial and test space
      typedef std::function<double(const VectorValues& u,
                                   const VectorValues& v)> VectorIntegrand;

      /// The same, depending on the quadrature point as well
      typedef std::function<double(const QuadPoint& p, const VectorValues& u,
                                   const VectorValues& v)> VectorPointIntegrand;

      /// Integrand with a vector valued trial and a scalar test function, e.g. -(div u, q)
      typedef std::function<double(const VectorValues& u,
                                   const PointValues& v)> VectorScalarIntegrand;

      /// Integrand with a scalar trial and a vector valued test function, e.g. -(p, div v)
      typedef std::function<double(const PointValues& u,
                                   const VectorValues& v)> ScalarVectorIntegrand;

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

      /// Adds the integral of the integrand over all cells
      void AddVolumeTerm(Integrand);

      /// Adds the integral of the integrand over all cells
      void AddVolumeTerm(PointIntegrand);

      /// Adds the integral of the integrand over all cells
      void AddVolumeTerm(VectorIntegrand);

      /// Adds the integral of the integrand over all cells
      void AddVolumeTerm(VectorPointIntegrand);

      /// Adds the integral of the integrand over all cells
      void AddVolumeTerm(VectorScalarIntegrand);

      /// Adds the integral of the integrand over all cells
      void AddVolumeTerm(ScalarVectorIntegrand);

      /**
       * Adds the integral of the integrand over the part of the boundary where the
       * indicator is true, over the whole boundary if it is omitted
       */
      void AddBoundaryTerm(Integrand, BoundaryIndicator = nullptr);

      /**
       * Adds the integral of the integrand over the part of the boundary where the
       * indicator is true, over the whole boundary if it is omitted
       */
      void AddBoundaryTerm(BoundaryIntegrand, BoundaryIndicator = nullptr);

      /**
       * Sets the values prescribed on the Dirichlet boundary. They enter the right hand side
       * as the lifting A_fd g_d, which DirichletRhs() returns after the assembly and which
       * has to be subtracted from the load vector. Without them the values are zero.
       */
      void SetDirichletValues(const DirichletValues&);

      /**
       * Assembles the finite element matrix.
       */
      void Assemble();

      /// Position of the matrix within a larger one, and whether it is transposed there
      struct Placement
      {
        size_t row, col;
        bool transposed;
      };

      /**
       * Assembles the matrix, or its transpose, into a larger one, starting at the given
       * row and column. Used for the blocks of a BlockSystem.
       */
      void Assemble(chemfem::linalg::SparseMatrixInserter&, size_t RowOffset, size_t ColOffset,
                    bool Transposed = false);

      /**
       * Assembles the matrix once and inserts it at several places of a larger one, e.g.
       * into a block and into its transpose.
       */
      void Assemble(chemfem::linalg::SparseMatrixInserter&, const std::vector<Placement>&);

      /**
       * Returns the matrix which corresponds to the bilinear form. Before calling this function
       * the assemble routine has to be invoked. Otherwise, an empty matrix is returned.
       */
      SparseMatrix& SystemMatrix();

      /**
       * The contribution of the prescribed Dirichlet values to the right hand side. Subtract
       * it from the load vector, it is zero without SetDirichletValues.
       */
      const Vector& DirichletRhs() const;

      const FESpace& GetTrialSpace() const;

      const FESpace& GetTestSpace() const;

    private:
      /// Adds the local matrix of a cell to the global one, only the free DOFs are kept
      void InsertLocalMatrix(chemfem::linalg::SparseMatrixInserter&,
                             const std::vector<Placement>&, size_t Cell,
                             const chemfem::linalg::DenseMatrix&);

      /// A boundary term, given by one of the two kinds of integrands
      struct BoundaryTerm
      {
        Integrand integrand;
        BoundaryIntegrand point_integrand;
        BoundaryIndicator part;
      };

      const FESpace& TrialSpace;
      const FESpace& TestSpace;

      SparseMatrix Matrix;
      Vector DirichletTerm;

      const DirichletValues* PrescribedValues = nullptr;

      std::vector<FEExpression> Terms;

      std::vector<Integrand> VolumeTerms;
      std::vector<PointIntegrand> PointVolumeTerms;
      std::vector<VectorIntegrand> VectorTerms;
      std::vector<VectorPointIntegrand> VectorPointTerms;
      std::vector<VectorScalarIntegrand> VectorScalarTerms;
      std::vector<ScalarVectorIntegrand> ScalarVectorTerms;
      std::vector<BoundaryTerm> BoundaryTerms;
    };

  };
};

#endif
