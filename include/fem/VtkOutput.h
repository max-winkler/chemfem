#ifndef _VTK_OUTPUT_H_
#define _VTK_OUTPUT_H_

#include <string>
#include <vector>

#include "fem/FEFunction.h"
#include "mesh/Mesh.h"

namespace chemfem{
  namespace fem{

    /**
     * Writes FE functions on a mesh into a legacy VTK file, e.g. for ParaView. Each field
     * is written with its values in the vertices. Where a function is discontinuous, as
     * for the Crouzeix-Raviart element, the values of the cells around a vertex are
     * averaged.
     *
     * The FE functions are only referenced, they have to live until Write() is called.
     */
    class VtkOutput
    {
    public:
      VtkOutput(const chemfem::mesh::Mesh&);

      /// Adds a scalar field under the given name, which must not contain spaces
      void AddScalar(const std::string&, const FEFunction&);

      /// Adds a vector field given by its two components
      void AddVector(const std::string&, const FEFunction&, const FEFunction&);

      /// Adds a vector field given by one FE function on a space with two components
      void AddVector(const std::string&, const FEFunction&);

      void Write(const std::string&) const;

    private:
      bool OnMesh(const FEFunction&) const;

      /// Values in the vertices of the given component of an FE function
      std::vector<double> VertexValues(const FEFunction&, int component = 0) const;

      struct ScalarField
      {
        std::string name;
        const FEFunction* function;
      };

      /**
       * A vector field, either given by two scalar FE functions or, if y is empty, by the
       * two components of x
       */
      struct VectorField
      {
        std::string name;
        const FEFunction *x, *y;
      };

      const chemfem::mesh::Mesh& mesh;

      std::vector<ScalarField> Scalars;
      std::vector<VectorField> Vectors;
    };

  };
};

#endif
