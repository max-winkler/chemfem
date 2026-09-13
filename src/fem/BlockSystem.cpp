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

      Values.assign(this->Spaces.size(), nullptr);
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

    void BlockSystem::SetDirichletValues(size_t i, const DirichletValues& g)
    {
      if(i >= Spaces.size() || &g.GetFESpace() != Spaces[i])
        {
          std::cerr << "Error: The Dirichlet values for the unknown " << i
                    << " have to belong to its space.\n";
          return;
        }

      Values[i] = &g;
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

    SparseMatrix& BlockSystem::AssembleMatrix()
    {
      const size_t n = Offset.back() + Constraints.size();

      Matrix = SparseMatrix(n, n);
      LU.reset();

      SparseMatrixInserter Ins(Matrix);

      // The values prescribed for an unknown belong to the forms of its column
      for(size_t b=0; b<Blocks.size(); ++b)
        {
          if(!Values[Blocks[b].col])
            continue;

          if(Blocks[b].transposed)
            {
              std::cerr << "Error: Prescribed values in the transposed block ("
                        << Blocks[b].row << "," << Blocks[b].col << ") are not supported.\n";
              continue;
            }

          Blocks[b].form->SetDirichletValues(*Values[Blocks[b].col]);
        }

      for(size_t a=0; a<Blocks.size(); ++a)
        for(size_t b=a+1; b<Blocks.size(); ++b)
          if(Blocks[a].form == Blocks[b].form && !Blocks[a].transposed && !Blocks[b].transposed
             && Blocks[a].col != Blocks[b].col
             && (Values[Blocks[a].col] || Values[Blocks[b].col]))
            std::cerr << "Error: The bilinear form of the blocks (" << Blocks[a].row << ","
                      << Blocks[a].col << ") and (" << Blocks[b].row << "," << Blocks[b].col
                      << ") carries prescribed values, but it can hold only one set. Use one "
                      << "form per column.\n";

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

      return Matrix;
    }

    Vector BlockSystem::AssembleRhs()
    {
      Vector Rhs(Offset.back() + Constraints.size());

      for(size_t b=0; b<RhsBlocks.size(); ++b)
        {
          LinearForm& form = *RhsBlocks[b].form;
          form.Assemble();

          const Vector& F = form.LoadVector();
          for(size_t k=0; k<F.size(); ++k)
            Rhs[Offset[RhsBlocks[b].row] + k] += F[k];
        }

      // The prescribed values move to the right hand side, A_fd g_d of every block
      for(size_t b=0; b<Blocks.size(); ++b)
        {
          if(Blocks[b].transposed || !Values[Blocks[b].col])
            continue;

          const Vector& lifting = Blocks[b].form->DirichletRhs();

          for(size_t k=0; k<lifting.size(); ++k)
            Rhs[Offset[Blocks[b].row] + k] -= lifting[k];
        }

      for(size_t d=0; d<FixedDofs.size(); ++d)
        Rhs[FixedDofs[d]] = 0.;

      return Rhs;
    }

    void BlockSystem::AddToRhs(Vector& Rhs, size_t i, const Vector& v) const
    {
      if(i >= Spaces.size() || v.size() != Spaces[i]->NrFreeDof())
        {
          std::cerr << "Error: The vector added to the right hand side of the block " << i
                    << " has the wrong length.\n";
          return;
        }

      for(size_t k=0; k<v.size(); ++k)
        Rhs[Offset[i] + k] += v[k];
    }

    Vector BlockSystem::FreeDof(size_t i, const Vector& X) const
    {
      Vector free(Spaces[i]->NrFreeDof());

      for(size_t k=0; k<free.size(); ++k)
        free[k] = X[Offset[i] + k];

      return free;
    }

    Vector BlockSystem::Solve(const Vector& Rhs)
    {
      if(!LU)
        LU = std::make_unique<DirectSolver>(Matrix, false);

      return LU->Solve(Rhs);
    }

    SparseMatrix& BlockSystem::SystemMatrix()
    {
      return Matrix;
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

      if(Values[i])
        u.CreateFunction(FreeDof, *Values[i]);
      else
        u.CreateFunction(FreeDof);

      return u;
    }

  }
}
