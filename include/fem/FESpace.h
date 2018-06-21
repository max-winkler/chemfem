#ifndef _FE_SPACE_H_
#define _FE_SPACE_H_

#include "fem/Element.h"
#include "mesh/Mesh.h"

using chemfem::mesh::Mesh;

namespace chemfem{
  namespace fem{   

    /**
     * This class represents a finite element space which is characterized by 
     * the mesh, the finite element type and the polynomial degree.
     */
    class FESpace
    {
      friend class BilinearForm;
      friend class LinearForm;
      
    public:
      /**
       * Initialize the finite element space by a reference to the mesh, the finite element 
       * type.
       */
      FESpace(Mesh&, Element&);

      /**
       * Returns the number of degrees of freedom.
       */
      size_t NrDof() const;

      /**
       * Returns the global index of a local degree of freedom. The first argument 
       * is the index of the cell, the second one the local index.
       */
      size_t GetGlobalIndex(size_t, size_t) const;

      /**
       * Returns the reference element.
       */
      const Element& RefElement() const;
      
    private:
      size_t *Dof;

      size_t nr_dof;
      
      int DofPerCell, DofPerEdge, IntDofPerCell;
      
      Element& refElement;
      Mesh& mesh;

      /**
       * Saves a list of the global indices where Dirichlet boundary conditions are imposed.
       */
      std::set<size_t> DirichletNodes;
      
      size_t *BdDof;

      /**
       * Creates the map between global and local degrees of freedom.
       */
      void CreateDofMap();

      /**
       * Creates a list of Dirichlet nodes according to the boundary conditions 
       * passed to FESpace in its contructor.
       */
      void CreateDirichletBcMap();
    };
    
  };
};

#endif
