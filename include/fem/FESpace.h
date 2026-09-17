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
      friend class FEFunction;

    public:
      /**
       * Initialize the finite element space by a reference to the mesh and the finite
       * element. Homogeneous Dirichlet conditions are imposed on the boundary edges whose
       * midpoint IsDirichlet accepts, every other boundary edge carries the natural
       * condition. Without an indicator nothing is fixed at all, so pass WholeBoundary to
       * put a homogeneous Dirichlet condition on the whole boundary.
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
       * True if the functions of this space are vector valued, either because the DOFs are
       * interleaved components as for a product element, or because the basis functions are
       * vectors already on the reference element as for Raviart-Thomas.
       */
      bool IsVectorValued() const;

      /**
       * True if all DOFs of the element sit in the interior of a cell, as for the DG
       * elements. No DOF of such a space touches the boundary, so essential conditions
       * cannot be imposed on it at all, and its functions are written as VTK cell data.
       */
      bool AllDofsInterior() const;

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
       * The element as a scalar one, or null if its shape functions are vectors. Resolved
       * once when the space is built, so the branch costs nothing in an assembly loop.
       */
      const ScalarElement* AsScalar() const;

      /// The element as a vector one, or null if its shape functions are scalars
      const VectorElement* AsVector() const;

      /**
       * Returns a reference to the finite element mesh.
       */
      const chemfem::mesh::Mesh& GetMesh() const;

      /**
       * Returns a pointer to the first element of the local Dof map for the i-th cell.
       */
      const size_t* GetLocalDofMap(size_t) const;

      /**
       * Sign of a local DOF, +1 or -1. Only the edge DOFs of an H(div) conforming element
       * carry -1, and only in the cell that traverses the edge against its own direction.
       */
      double LocalSign(size_t cell, size_t local) const;

      /// True if the cell runs through its local edge against that edge's own direction
      bool IsEdgeReversed(size_t cell, int local_edge) const;

      /**
       * The frame the derivative DOFs of a vertex are taken along, see NodeFrame. Asked for
       * by the cell and the local vertex, because the DOF transformation works per cell.
       */
      NodeFrame VertexFrame(size_t cell, int local_vertex) const;

      /**
       * Interpolates a smooth function into the finite element space
       */
      FEFunction Interpolate(ScalarFunction);

      /**
       * The same for a vector valued space: each DOF takes the component of the value that
       * belongs to it, e.g. to lay an initial velocity into the space.
       */
      FEFunction Interpolate(VectorFunction);

    private:
      /**
       * Expands a vector of the free DOFs to one over all of them, with zero or with the
       * prescribed values at the Dirichlet DOFs. This is internal mechanics of
       * FEFunction::CreateFunction, which is how an FE function is built from a solution.
       */
      Vector IncorporateBC(const Vector&) const;
      Vector IncorporateBC(const Vector&, const DirichletValues&) const;

      Element& refElement;

      const ScalarElement* scalar;
      const VectorElement* vector;

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
