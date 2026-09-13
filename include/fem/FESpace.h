#ifndef _FE_SPACE_H_
#define _FE_SPACE_H_

#include "fem/DofManager.h"
#include "fem/Element.h"
#include "fem/FEExpression.h"
#include "mesh/Mesh.h"

namespace chemfem{
  namespace fem{

    // Forward declarations
    class FEFunction;
    class DirichletValues;

    /**
     * This class represents a finite element space which is characterized by
     * the mesh, the finite element type and the polynomial degree.
     */
    class FESpace
    {
      friend class BilinearForm;
      friend class LinearForm;
      friend class DirichletValues;

    public:
      /**
       * Initialize the finite element space by a reference to the mesh and the finite
       * element. Homogeneous Dirichlet conditions are imposed on the boundary edges whose
       * midpoint IsDirichlet accepts, on the whole boundary if it is omitted. All other
       * boundary edges are Neumann edges.
       */
      FESpace(chemfem::mesh::Mesh&, Element&, BoundaryIndicator IsDirichlet = nullptr);

      /**
       * Returns the number of degrees of freedom.
       */
      size_t NrDof() const;

      /**
       * Returns the number of free DOFs (DOFs that are not fixed by essential boundary conditions).
       */
      size_t NrFreeDof() const;

      /**
       * Returns the number of local degrees of freedom
       */
      size_t NrLocalDof() const;

      /// Number of components of the FE functions, 1 for a scalar space
      size_t NrComponents() const;

      /**
       * Returns the global index of a local degree of freedom. The first argument
       * is the index of the cell, the second one the local index.
       */
      size_t GetGlobalIndex(size_t, size_t) const;

      /**
       * Returns the reference element.
       */
      const Element& RefElement() const;

      /**
       * For a given vector representing the values at the free DOFs, this function returns
       * a vector where also Dirichlet DOFs are incorporated.
       */
      Vector IncorporateBC(const Vector&) const;

      /**
       * The same, with the prescribed values at the Dirichlet DOFs instead of zero.
       */
      Vector IncorporateBC(const Vector&, const DirichletValues&) const;

      /**
       * Returns a reference to the finite element mesh.
       */
      const chemfem::mesh::Mesh& GetMesh() const;

      /**
       * Returns a pointer to the first element of the local Dof map for the i-th cell.
       */
      const size_t* GetLocalDofMap(size_t) const;

      /**
       * Interpolates a smooth function into the finite element space
       */
      FEFunction Interpolate(ScalarFunction);

    private:
      Element& refElement;

      chemfem::mesh::Mesh& mesh;

      DofManager Dofs;

      /**
       * Write some information on the finite element space to a stream
       */
      friend std::ostream& operator<<(std::ostream&, const FESpace&);
    };

  };
};

#endif
