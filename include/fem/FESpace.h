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
       * Returns the number of free DOFs (DOFs that are not fixed by essential boundary conditions).
       */
      size_t NrFreeDof() const;
      
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
       * a vector where also Dirichlet DOFs are invorporated.
       */
      Vector IncorporateBC(const Vector&) const; 

      /**
       * Returns a reference to the finite element mesh.
       */
      const Mesh& GetMesh() const;

      /**
       * Returns a pointer to the first element of the local Dof map for the i-th cell.
       */
      const size_t* GetLocalDofMap(size_t) const;
      
    private:
      /**
       * Contains the mappings between local and global DOFs for each element
       */
      size_t *DofMap;
      
      /**
       * List of booleans. If the i-th index is true, the degree of freedom is
       * free, this means, no essential boundary conditions are present there.
       */
      bool *DofType;

      /**
       * This array is used to map all free and active degrees of freedom to a 
       * new numbering.
       */
      size_t *DofIndex;

      size_t nr_dof;
      size_t nr_free_dof;
      
      int DofPerCell, DofPerEdge, IntDofPerCell;
      
      Element& refElement;

      Mesh& mesh;

      /**
       * Saves a list of the global indices where Dirichlet boundary conditions are imposed.
       */
      std::set<size_t> DirichletNodes;
      
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
