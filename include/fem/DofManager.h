#ifndef _DOF_MANAGER_H_
#define _DOF_MANAGER_H_

#include <vector>

#include "fem/Element.h"
#include "fem/FEExpression.h"
#include "mesh/Mesh.h"

namespace chemfem{
  namespace fem{

    /**
     * Numbers the degrees of freedom of a finite element space and separates the free
     * ones from those fixed by Dirichlet conditions.
     *
     * The element states how many DOFs sit on each vertex, on each edge and in the
     * interior of a cell. Globally the vertex DOFs come first, then the edge DOFs, then
     * the interior DOFs. A cell orders its local DOFs the same way, the edge DOFs by
     * the local edge k, which runs from vertex k to vertex k+1.
     */
    class DofManager
    {
    public:
      /**
       * Numbers the DOFs of the element on the mesh. Dirichlet conditions are imposed
       * on the boundary edges whose midpoint IsDirichlet accepts, on the whole boundary
       * if it is empty.
       */
      DofManager(const chemfem::mesh::Mesh&, const Element&, const BoundaryIndicator&);

      /// Number of all DOFs
      size_t NrDof() const;

      /// Number of the DOFs that are not fixed by Dirichlet conditions
      size_t NrFreeDof() const;

      /// Number of DOFs per cell
      size_t NrLocalDof() const;

      /// Global index of a local DOF, given by the cell and its local index
      size_t GlobalIndex(size_t, size_t) const;

      /// Global indices of all local DOFs of a cell
      const size_t* LocalDofMap(size_t) const;

      /**
       * Sign of a local DOF, +1 or -1. An edge DOF of an H(div) conforming element measures
       * the flux through a globally fixed normal of the edge, while the reference basis
       * function refers to the outward normal of its own cell. The two neighbours of an edge
       * see opposite normals, so the cell that traverses the edge against its own direction
       * gets -1. Everything else is +1.
       */
      double LocalSign(size_t cell, size_t local) const;

      /// False if the DOF is fixed by a Dirichlet condition
      bool IsFree(size_t) const;

      /// Index of a free DOF among the free ones, of a Dirichlet DOF among the Dirichlet ones
      size_t ReducedIndex(size_t) const;

      /// True for the boundary edges with a Dirichlet condition
      bool IsDirichletEdge(size_t) const;

    private:
      void CreateDofMap(const chemfem::mesh::Mesh&, const Element&);

      void MarkDirichletDofs(const chemfem::mesh::Mesh&, const Element&,
                             const BoundaryIndicator&);

      size_t nr_local_dof, nr_dof, nr_free_dof;

      /// Global indices of the local DOFs, NrLocalDof() consecutive entries per cell
      std::vector<size_t> DofMap;

      /// Signs of the local DOFs, laid out like DofMap. Empty unless the element needs them.
      std::vector<double> Sign;

      std::vector<bool> Free;
      std::vector<size_t> Reduced;
      std::vector<bool> DirichletEdge;
    };

  };
};

#endif
