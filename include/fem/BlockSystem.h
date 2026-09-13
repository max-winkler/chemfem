#ifndef _BLOCK_SYSTEM_H_
#define _BLOCK_SYSTEM_H_

#include <functional>
#include <initializer_list>
#include <memory>
#include <vector>

#include "fem/BilinearForm.h"
#include "fem/LinearForm.h"
#include "fem/FEFunction.h"
#include "linalg/DirectSolver.h"
#include "linalg/SparseMatrix.h"
#include "linalg/Vector.h"

namespace chemfem{
  namespace fem{

    /**
     * Linear system for several coupled unknowns, each one in its own FE space, e.g. the
     * two velocity components and the pressure of the Stokes equations:
     *
     * \code
     *   BlockSystem S({V, V, Q});
     *   S.AddBlock(0, 0, A);   S.AddBlock(1, 1, A);
     *   S.AddBlock(2, 0, Bx);  S.AddTransposedBlock(0, 2, Bx);
     *   S.AddRhs(0, Fx);
     *   S.FixDof(2);
     *   S.AssembleMatrix();
     *
     *   FEFunction Ux = S.Extract(0, S.Solve(S.AssembleRhs()));
     * \endcode
     *
     * The unknowns are numbered block by block. The block (i,j) couples the test functions
     * of the space i with the trial functions of the space j and is given by bilinear
     * forms with exactly these spaces. A form used in several blocks is assembled only once.
     */
    class BlockSystem
    {
    public:
      BlockSystem(std::initializer_list<std::reference_wrapper<const FESpace> >);

      /// Adds a bilinear form to the block (i,j), several forms on one block add up
      void AddBlock(size_t i, size_t j, BilinearForm&);

      /**
       * Adds the transpose of a bilinear form to the block (i,j). The form acts on the test
       * space j and the trial space i. For the Stokes equations, the divergence -(div u, q)
       * gives the pressure term -(p, div v) this way.
       */
      void AddTransposedBlock(size_t i, size_t j, BilinearForm&);

      /// Adds a linear form to the right hand side of the block i
      void AddRhs(size_t i, LinearForm&);

      /**
       * Prescribes the Dirichlet values of the unknown i. They enter the right hand side as
       * the lifting of the blocks in the column i, and Extract puts them into the solution.
       * A form that is used in several columns cannot carry values, because it holds only
       * one set of them. Use one form per column in that case.
       */
      void SetDirichletValues(size_t i, const DirichletValues&);

      /**
       * Fixes one degree of freedom of the unknown i to zero, e.g. to determine the constant
       * of a pressure. This keeps the matrix sparse, unlike AddMeanValueConstraint, and the
       * solution is normalized afterwards with FEFunction::SubtractMean.
       */
      void FixDof(size_t i, size_t dof = 0);

      /**
       * Fixes the constant in the unknown i by requiring that its integral vanishes. The
       * constraint is imposed with a Lagrange multiplier, which adds one row and one column
       * to the system. It couples all DOFs of that unknown, which makes the factorization of
       * a direct solver considerably more expensive than FixDof.
       */
      void AddMeanValueConstraint(size_t i);

      /// Assembles the matrix from its blocks and constraints
      SparseMatrix& AssembleMatrix();

      /**
       * Assembles the right hand side from the linear forms and the Dirichlet lifting. It does
       * not depend on the matrix, so with data that does not change in time it can be computed
       * once outside a time loop.
       */
      Vector AssembleRhs();

      /**
       * Adds a vector to the block i of a right hand side. Used for contributions that are not
       * an integral over the mesh, e.g. the mass matrix times the solution of the previous
       * time step.
       */
      void AddToRhs(Vector&, size_t i, const Vector&) const;

      /// The free degrees of freedom of the unknown i, taken from the solution of the system
      Vector FreeDof(size_t i, const Vector&) const;

      /**
       * Solves the system for the given right hand side with UMFPACK, without its iterative
       * refinement. The factorization is computed on the first call and kept until the matrix
       * is assembled again, so further right hand sides only cost a forward and a backward
       * substitution. For an ill conditioned system, build a DirectSolver from SystemMatrix().
       */
      Vector Solve(const Vector&);

      SparseMatrix& SystemMatrix();

      /// Number of FE unknowns, without the Lagrange multipliers of the constraints
      size_t NrDof() const;

      /// The FE function of the unknown i, taken from the solution of the whole system
      FEFunction Extract(size_t i, const Vector&) const;

    private:
      struct Block
      {
        size_t row, col;
        BilinearForm* form;
        bool transposed;
      };

      struct RhsBlock
      {
        size_t row;
        LinearForm* form;
      };

      std::vector<const FESpace*> Spaces;

      /// First row of each block, the last entry is the number of all FE unknowns
      std::vector<size_t> Offset;

      std::vector<Block> Blocks;
      std::vector<RhsBlock> RhsBlocks;

      /// Unknowns with a mean value constraint, one Lagrange multiplier each
      std::vector<size_t> Constraints;

      /// Indices of the unknowns fixed by FixDof
      std::vector<size_t> FixedDofs;

      /// The prescribed values of each unknown, empty where there are none
      std::vector<const DirichletValues*> Values;

      SparseMatrix Matrix;

      std::unique_ptr<chemfem::linalg::DirectSolver> LU;
    };

  };
};

#endif
