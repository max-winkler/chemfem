#include "fem/VtkOutput.h"

#include "fem/DofTransform.h"

#include "Verbosity.h"

#include <fstream>
#include <iomanip>
#include <iostream>

namespace chemfem{
  namespace fem{

    using chemfem::mesh::Mesh;

    VtkOutput::VtkOutput(const Mesh& mesh) : mesh(mesh) {}

    bool VtkOutput::OnMesh(const FEFunction& u) const
    {
      if(&u.GetFESpace().GetMesh() == &mesh)
        return true;

      std::cerr << "Error: The FE function does not live on the mesh of the VTK output.\n";
      return false;
    }

    void VtkOutput::AddScalar(const std::string& name, const FEFunction& u)
    {
      if(OnMesh(u))
        Scalars.push_back(ScalarField{name, &u});
    }

    void VtkOutput::AddVector(const std::string& name, const FEFunction& ux,
                              const FEFunction& uy)
    {
      if(OnMesh(ux) && OnMesh(uy))
        Vectors.push_back(VectorField{name, &ux, &uy});
    }

    void VtkOutput::AddVector(const std::string& name, const FEFunction& u)
    {
      // Either two interleaved components, as for a product element, or a Piola mapped
      // element, whose basis functions are vectors in themselves
      const FESpace& Space = u.GetFESpace();

      if(Space.NrComponents() != 2 && !Space.AsVector())
        {
          std::cerr << "Error: The vector field " << name << " needs an FE function with two "
                    << "components.\n";
          return;
        }

      if(OnMesh(u))
        Vectors.push_back(VectorField{name, &u, nullptr});
    }

    std::vector<double> VtkOutput::VertexValues(const FEFunction& u, int component) const
    {
      const FESpace& Space = u.GetFESpace();
      const double Vertex[3][2] = {{0., 0.}, {1., 0.}, {0., 1.}};

      std::vector<double> sum(mesh.NrNodes(), 0.);
      std::vector<int> count(mesh.NrNodes(), 0);

      // A Piola mapped basis function contributes to both components, so the filter on the
      // component does not apply to it and its value needs the mapping of its cell
      const VectorElement* Vec = Space.AsVector();
      const bool Piola = Vec != nullptr;

      std::vector<double> Coeff;

      for(size_t c=0; c<mesh.NrCells(); ++c)
      {
        GatherLocalCoefficients(Space, u.Coefficients(), c, Coeff);

        for(int v=0; v<3; ++v)
          {
            double value = 0.;
            for(size_t k=0; k<Space.NrLocalDof(); ++k)
              {
                if(!Piola && Space.RefElement().Component(k) != component)
                  continue;

                const double basis = Piola
                  ? MapFromReference(ReferenceVector{Vec->Value(k, Vertex[v][0], Vertex[v][1]),
                                                     Vec->Gradient(k, Vertex[v][0], Vertex[v][1])},
                                     mesh.Jacobian(c), mesh.Determinant(c),
                                     Space.LocalSign(c, k)).value[component]
                  : Space.AsScalar()->Value(k, Vertex[v][0], Vertex[v][1]);

                value += Coeff[k] * basis;
              }

            const size_t node = mesh.Cells[c].LocNode[v];
            sum[node] += value;
            ++count[node];
          }
      }

      for(size_t n=0; n<sum.size(); ++n)
        if(count[n] > 0)
          sum[n] /= count[n];

      return sum;
    }

    bool VtkOutput::IsCellData(const FEFunction& u) const
    {
      return u.GetFESpace().AllDofsInterior();
    }

    std::vector<double> VtkOutput::CellValues(const FEFunction& u, int component) const
    {
      const FESpace& Space = u.GetFESpace();

      std::vector<double> values(mesh.NrCells(), 0.);
      std::vector<double> Coeff;

      for(size_t c=0; c<mesh.NrCells(); ++c)
        {
          GatherLocalCoefficients(Space, u.Coefficients(), c, Coeff);

          double value = 0.;

          for(size_t k=0; k<Space.NrLocalDof(); ++k)
            {
              if(Space.RefElement().Component(k) != component)
                continue;

              value += Coeff[k]
                * Space.AsScalar()->Value(k, 1./3, 1./3);
            }

          values[c] = value;
        }

      return values;
    }

    void VtkOutput::Write(const std::string& filename) const
    {
      std::ofstream ofs(filename);

      if(!ofs)
        {
          std::cerr << "Error: Cannot open " << filename << " for writing.\n";
          return;
        }

      ofs << std::setprecision(12);

      ofs << "# vtk DataFile Version 3.0\n"
          << "chemfem\n"
          << "ASCII\n"
          << "DATASET UNSTRUCTURED_GRID\n";

      ofs << "POINTS " << mesh.NrNodes() << " double\n";
      for(size_t n=0; n<mesh.NrNodes(); ++n)
        ofs << mesh.Nodes[n].getX() << " " << mesh.Nodes[n].getY() << " 0\n";

      ofs << "CELLS " << mesh.NrCells() << " " << 4*mesh.NrCells() << "\n";
      for(size_t c=0; c<mesh.NrCells(); ++c)
        ofs << "3 " << mesh.Cells[c].LocNode[0] << " " << mesh.Cells[c].LocNode[1]
            << " " << mesh.Cells[c].LocNode[2] << "\n";

      // 5 is the VTK cell type of a triangle
      ofs << "CELL_TYPES " << mesh.NrCells() << "\n";
      for(size_t c=0; c<mesh.NrCells(); ++c)
        ofs << "5\n";

      if(Scalars.empty() && Vectors.empty())
        {
          if(Verbose())
            std::cout << "Wrote " << filename << std::endl;
          return;
        }

      std::vector<size_t> PointScalars, CellScalars, PointVectors, CellVectors;

      for(size_t s=0; s<Scalars.size(); ++s)
        (IsCellData(*Scalars[s].function) ? CellScalars : PointScalars).push_back(s);

      for(size_t v=0; v<Vectors.size(); ++v)
        (IsCellData(*Vectors[v].x) ? CellVectors : PointVectors).push_back(v);

      if(!PointScalars.empty() || !PointVectors.empty())
        {
          ofs << "POINT_DATA " << mesh.NrNodes() << "\n";

          for(size_t i=0; i<PointScalars.size(); ++i)
            {
              const ScalarField& field = Scalars[PointScalars[i]];
              const std::vector<double> values = VertexValues(*field.function);

              ofs << "SCALARS " << field.name << " double 1\n"
                  << "LOOKUP_TABLE default\n";
              for(size_t n=0; n<values.size(); ++n)
                ofs << values[n] << "\n";
            }

          for(size_t i=0; i<PointVectors.size(); ++i)
            {
              const VectorField& field = Vectors[PointVectors[i]];
              const std::vector<double> x = VertexValues(*field.x, 0);
              const std::vector<double> y = field.y ? VertexValues(*field.y, 0)
                                                    : VertexValues(*field.x, 1);

              ofs << "VECTORS " << field.name << " double\n";
              for(size_t n=0; n<x.size(); ++n)
                ofs << x[n] << " " << y[n] << " 0\n";
            }
        }

      if(!CellScalars.empty() || !CellVectors.empty())
        {
          ofs << "CELL_DATA " << mesh.NrCells() << "\n";

          for(size_t i=0; i<CellScalars.size(); ++i)
            {
              const ScalarField& field = Scalars[CellScalars[i]];
              const std::vector<double> values = CellValues(*field.function);

              ofs << "SCALARS " << field.name << " double 1\n"
                  << "LOOKUP_TABLE default\n";
              for(size_t c=0; c<values.size(); ++c)
                ofs << values[c] << "\n";
            }

          for(size_t i=0; i<CellVectors.size(); ++i)
            {
              const VectorField& field = Vectors[CellVectors[i]];
              const std::vector<double> x = CellValues(*field.x, 0);
              const std::vector<double> y = field.y ? CellValues(*field.y, 0)
                                                    : CellValues(*field.x, 1);

              ofs << "VECTORS " << field.name << " double\n";
              for(size_t c=0; c<x.size(); ++c)
                ofs << x[c] << " " << y[c] << " 0\n";
            }
        }

      if(Verbose())
        std::cout << "Wrote " << filename << std::endl;
    }

  }
}
