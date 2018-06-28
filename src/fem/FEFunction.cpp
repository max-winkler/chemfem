#include "fem/FEFunction.h"

namespace chemfem{
  namespace fem{
    
    FEFunction::FEFunction(const FEFunction& other)
      : Space(other.Space), Data(other.Data) {}

    FEFunction::FEFunction(const FESpace& Space)
      : Space(&Space), Data(Space.NrDof()) {}

    const FESpace& FEFunction::GetFESpace() const
    {
      return *Space;
    }

    const double& FEFunction::operator[](size_t k) const
    {
      return Data[k];
    }

    void FEFunction::CreateFunction(const Vector& FreeDof)
    {
      Data = Space->IncorporateBC(FreeDof);
    }

    void FEFunction::SetCoefficients(const Vector& Data)
    {
      this->Data = Data;
    }

    FEFunction FESpace::Interpolate(double (*u)(double, double))
    {
      FEFunction Function(*this);
      Vector Vec(nr_dof);
      
      // TODO Set functions
      Function.SetCoefficients(Vec);
      return Function;
    }

    void FEFunction::WriteVtk(const std::string& filename) const
    {
      Space->GetMesh().WriteVtk(filename, Data);
    }
  }
}
