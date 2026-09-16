#include "fem/LinearSystem.h"

namespace chemfem{
  namespace fem{

    LinearSystem::LinearSystem(const FESpace& Space) : System({Space}) {}

    void LinearSystem::AddLhs(BilinearForm& a)
    {
      System.AddBlock(0, 0, a);
    }

    void LinearSystem::AddRhs(LinearForm& l)
    {
      System.AddRhs(0, l);
    }

    void LinearSystem::SetDirichletValues(const DirichletValues& g)
    {
      System.SetDirichletValues(0, g);
    }

    void LinearSystem::FixDof(size_t dof)
    {
      System.FixDof(0, dof);
    }

    void LinearSystem::AddMeanValueConstraint()
    {
      System.AddMeanValueConstraint(0);
    }

    SparseMatrix& LinearSystem::AssembleMatrix()
    {
      return System.AssembleMatrix();
    }

    Vector LinearSystem::AssembleRhs()
    {
      return System.AssembleRhs();
    }

    void LinearSystem::AddToRhs(Vector& Rhs, const Vector& contribution) const
    {
      System.AddToRhs(Rhs, 0, contribution);
    }

    Vector LinearSystem::FreeDof(const Vector& Solution) const
    {
      return System.FreeDof(0, Solution);
    }

    Vector LinearSystem::Solve(const Vector& Rhs)
    {
      return System.Solve(Rhs);
    }

    SparseMatrix& LinearSystem::SystemMatrix()
    {
      return System.SystemMatrix();
    }

    size_t LinearSystem::NrDof() const
    {
      return System.NrDof();
    }

    FEFunction LinearSystem::Extract(const Vector& Solution) const
    {
      return System.Extract(0, Solution);
    }

  }
}
