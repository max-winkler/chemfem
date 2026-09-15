#include "fem/VtkOutput.h"

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
      if(u.GetFESpace().NrComponents() != 2)
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

      for(size_t c=0; c<mesh.NrCells(); ++c)
        for(int v=0; v<3; ++v)
          {
            double value = 0.;
            for(size_t k=0; k<Space.NrLocalDof(); ++k)
              {
                if(Space.RefElement().Component(k) != component)
                  continue;

                value += u[Space.GetGlobalIndex(c, k)]
                  * Space.RefElement().Value(k, Vertex[v][0], Vertex[v][1]);
              }

            const size_t node = mesh.Cells[c].LocNode[v];
            sum[node] += value;
            ++count[node];
          }

      for(size_t n=0; n<sum.size(); ++n)
        if(count[n] > 0)
          sum[n] /= count[n];

      return sum;
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

      ofs << "POINT_DATA " << mesh.NrNodes() << "\n";

      for(size_t s=0; s<Scalars.size(); ++s)
        {
          const std::vector<double> values = VertexValues(*Scalars[s].function);

          ofs << "SCALARS " << Scalars[s].name << " double 1\n"
              << "LOOKUP_TABLE default\n";
          for(size_t n=0; n<values.size(); ++n)
            ofs << values[n] << "\n";
        }

      for(size_t v=0; v<Vectors.size(); ++v)
        {
          const std::vector<double> x = VertexValues(*Vectors[v].x, 0);
          const std::vector<double> y = Vectors[v].y ? VertexValues(*Vectors[v].y, 0)
                                                     : VertexValues(*Vectors[v].x, 1);

          ofs << "VECTORS " << Vectors[v].name << " double\n";
          for(size_t n=0; n<x.size(); ++n)
            ofs << x[n] << " " << y[n] << " 0\n";
        }

      if(Verbose())
        std::cout << "Wrote " << filename << std::endl;
    }

  }
}
