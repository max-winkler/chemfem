#ifndef _BLOCK_SYSTEM_H_
#define _BLOCK_SYSTEM_H_

#include <vector>

#include "fem/BilinearForm.h"
#include "fem/LinearForm.h"
#include "fem/FEFunction.h"
#include "linalg/SparseMatrix.h"
#include "linalg/Vector.h"

namespace chemfem{
  namespace fem{

    /**
     * Linear system for several coupled unknowns, each one in its own FE space, e.g. the
     * two velocity components and the pressure of the Stokes equations.
     *
     * The unknowns are numbered block by block. The block (i,j) couples the test functions
     * of the space i with the trial functions of the space j and is given by bilinear
     * forms with exactly these spaces.
     */
    class BlockSystem
    {
    public:
      BlockSystem(const std::vector<const FESpace*>&);

      /// Adds a bilinear form to the block (i,j), several forms on one block add up
      void AddBlock(size_t i, size_t j, BilinearForm&);

      /// Adds a linear form to the right hand side of the block i
      void AddRhs(size_t i, LinearForm&);

      /**
       * Fixes the constant in the unknown i, e.g. the pressure, by requiring that its
       * integral vanishes. The constraint is imposed with a Lagrange multiplier, which
       * adds one row and one column to the system.
       */
      void AddMeanValueConstraint(size_t i);

      /// Assembles the matrix and the right hand side
      void Assemble();

      /// Assembles only the matrix, with its blocks and constraints
      void AssembleMatrix();

      /**
       * Assembles only the right hand side. In a time stepping scheme with a constant step
       * size the matrix stays the same, and only this has to be redone in each step.
       */
      void AssembleRhs();

      SparseMatrix& SystemMatrix();

      Vector& Rhs();

      /// The FE function of the unknown i, taken from the solution of the whole system
      FEFunction Extract(size_t i, const Vector&) const;

    private:
      struct Block
      {
        size_t row, col;
        BilinearForm* form;
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

      SparseMatrix Matrix;
      Vector RhsVector;
    };

  };
};

#endif
