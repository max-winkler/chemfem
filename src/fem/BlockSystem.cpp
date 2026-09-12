#include "fem/BlockSystem.h"

#include <iostream>

#include "linalg/SparseMatrixInserter.h"

using chemfem::linalg::Coordinate;
using chemfem::linalg::DirectSolver;
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

    BlockSystem::BlockSystem(std::initializer_list<std::reference_wrapper<const FESpace> > Spaces)
      : Offset(Spaces.size()+1, 0), Matrix(0, 0)
    {
      for(std::reference_wrapper<const FESpace> space : Spaces)
        this->Spaces.push_back(&space.get());

      for(size_t i=0; i<this->Spaces.size(); ++i)
        Offset[i+1] = Offset[i] + this->Spaces[i]->NrFreeDof();
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

      Blocks.push_back(Block{i, j, &form, false});
    }

    void BlockSystem::AddTransposedBlock(size_t i, size_t j, BilinearForm& form)
    {
      if(i >= Spaces.size() || j >= Spaces.size()
         || &form.GetTestSpace() != Spaces[j] || &form.GetTrialSpace() != Spaces[i])
        {
          std::cerr << "Error: The bilinear form for the transposed block (" << i << "," << j
                    << ") needs the test space " << j << " and the trial space " << i
                    << ".\n";
          return;
        }

      Blocks.push_back(Block{i, j, &form, true});
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

    void BlockSystem::FixDof(size_t i, size_t dof)
    {
      if(i >= Spaces.size() || dof >= Spaces[i]->NrFreeDof())
        {
          std::cerr << "Error: The unknown " << i << " has no free DOF " << dof << ".\n";
          return;
        }

      FixedDofs.push_back(Offset[i] + dof);
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
      LU.reset();

      SparseMatrixInserter Ins(Matrix);

      // Every form is assembled once and inserted into all blocks it appears in
      std::vector<BilinearForm*> Forms;
      std::vector<std::vector<BilinearForm::Placement> > Places;

      for(size_t b=0; b<Blocks.size(); ++b)
        {
          size_t f = 0;
          while(f < Forms.size() && Forms[f] != Blocks[b].form)
            ++f;

          if(f == Forms.size())
            {
              Forms.push_back(Blocks[b].form);
              Places.push_back(std::vector<BilinearForm::Placement>());
            }

          Places[f].push_back(BilinearForm::Placement{Offset[Blocks[b].row],
                                                      Offset[Blocks[b].col],
                                                      Blocks[b].transposed});
        }

      for(size_t f=0; f<Forms.size(); ++f)
        Forms[f]->Assemble(Ins, Places[f]);

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

      // The diagonal entry of a fixed DOF has to be in the sparsity pattern
      for(size_t d=0; d<FixedDofs.size(); ++d)
        Ins.Insert(FixedDofs[d], FixedDofs[d], 1.);

      Ins.Build();

      for(size_t d=0; d<FixedDofs.size(); ++d)
        Matrix.EliminateRowAndColumn(FixedDofs[d]);
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

      for(size_t d=0; d<FixedDofs.size(); ++d)
        RhsVector[FixedDofs[d]] = 0.;
    }

    Vector BlockSystem::Solve(bool IterativeRefinement)
    {
      if(!LU)
        LU = std::make_unique<DirectSolver>(Matrix, IterativeRefinement);

      return LU->Solve(RhsVector);
    }

    SparseMatrix& BlockSystem::SystemMatrix()
    {
      return Matrix;
    }

    Vector& BlockSystem::Rhs()
    {
      return RhsVector;
    }

    size_t BlockSystem::NrDof() const
    {
      return Offset.back();
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
