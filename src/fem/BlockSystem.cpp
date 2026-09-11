#include "fem/BlockSystem.h"

#include <iostream>

#include "linalg/SparseMatrixInserter.h"

using chemfem::linalg::Coordinate;
using chemfem::linalg::SparseMatrix;
using chemfem::linalg::SparseMatrixInserter;
using chemfem::linalg::Vector;

namespace chemfem{
  namespace fem{

    namespace {

      double One(const Coordinate&)
      {
        return 1.;
      }
    }

    BlockSystem::BlockSystem(const std::vector<const FESpace*>& Spaces)
      : Spaces(Spaces), Offset(Spaces.size()+1, 0), Matrix(0, 0)
    {
      for(size_t i=0; i<Spaces.size(); ++i)
        Offset[i+1] = Offset[i] + Spaces[i]->NrFreeDof();
    }

    void BlockSystem::AddBlock(size_t i, size_t j, BilinearForm& form)
    {
      if(i >= Spaces.size() || j >= Spaces.size()
         || &form.GetTestSpace() != Spaces[i] || &form.GetTrialSpace() != Spaces[j])
        {
          std::cerr << "Error: The bilinear form for the block (" << i << "," << j
                    << ") needs the test space " << i << " and the trial space " << j
                    << ".\n";
          return;
        }

      Blocks.push_back(Block{i, j, &form});
    }

    void BlockSystem::AddRhs(size_t i, LinearForm& form)
    {
      if(i >= Spaces.size() || &form.GetTestSpace() != Spaces[i])
        {
          std::cerr << "Error: The linear form for the block " << i
                    << " needs the test space " << i << ".\n";
          return;
        }

      RhsBlocks.push_back(RhsBlock{i, &form});
    }

    void BlockSystem::AddMeanValueConstraint(size_t i)
    {
      Constraints.push_back(i);
    }

    void BlockSystem::Assemble()
    {
      AssembleMatrix();
      AssembleRhs();
    }

    void BlockSystem::AssembleMatrix()
    {
      const size_t n = Offset.back() + Constraints.size();

      Matrix = SparseMatrix(n, n);

      SparseMatrixInserter Ins(Matrix);

      for(size_t b=0; b<Blocks.size(); ++b)
        Blocks[b].form->Assemble(Ins, Offset[Blocks[b].row], Offset[Blocks[b].col]);

      // The multiplier couples with the integrals of the basis functions of its unknown
      for(size_t c=0; c<Constraints.size(); ++c)
        {
          const size_t i = Constraints[c];
          const size_t row = Offset.back() + c;

          LinearForm Mean(*Spaces[i]);
          Mean.AddVolumeForce(One);
          Mean.Assemble();

          const Vector& m = Mean.LoadVector();
          for(size_t k=0; k<m.size(); ++k)
            {
              Ins.Insert(row, Offset[i] + k, m[k]);
              Ins.Insert(Offset[i] + k, row, m[k]);
            }
        }

      Ins.Build();
    }

    void BlockSystem::AssembleRhs()
    {
      RhsVector = Vector(Offset.back() + Constraints.size());

      for(size_t b=0; b<RhsBlocks.size(); ++b)
        {
          LinearForm& form = *RhsBlocks[b].form;
          form.Assemble();

          const Vector& F = form.LoadVector();
          for(size_t k=0; k<F.size(); ++k)
            RhsVector[Offset[RhsBlocks[b].row] + k] += F[k];
        }
    }

    SparseMatrix& BlockSystem::SystemMatrix()
    {
      return Matrix;
    }

    Vector& BlockSystem::Rhs()
    {
      return RhsVector;
    }

    FEFunction BlockSystem::Extract(size_t i, const Vector& X) const
    {
      Vector FreeDof(Spaces[i]->NrFreeDof());

      for(size_t k=0; k<FreeDof.size(); ++k)
        FreeDof[k] = X[Offset[i] + k];

      FEFunction u(*Spaces[i]);
      u.CreateFunction(FreeDof);

      return u;
    }

  }
}
