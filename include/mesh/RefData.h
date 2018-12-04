#ifndef _REF_DATA_H_
#define _REF_DATA_H_

namespace chemfem{
  namespace mesh{

    enum RefType{NO_REF, REGULAR_REF};
    
    /**
     * This is an abstract class which is used for mesh refinement.
     * In all derived classes the information needed to refine a 
     * single cell are stored.
     */
    class RefData
    {
    public:

      /// Constructor which initailized the number of new cells, nodes and edges
      RefData(int, int, int);
      
      /// Returns the number of nodes
      int GetNrNodes();
      /// Returns the number of edges
      int GetNrEdges();
      /// Returns the number of cells
      int GetNrCells();
      /// Returns a pointer to the coordinates of the i-th new node
      double* GetNodeCoords(int);
      /// Returns a pointer to the endpoints of the edges (dimension is 2)
      int* GetEdge(int);
      /// Returns a pointer to the vertices of the new cells (dimension is 3)
      int* GetCell(int);
      
    protected:
      /// The number of new cells (after refinement)
      int NrNewCells;
      /// The indices of the vertices of new cells
      int* NewCells;            
      /// The number of new nodes (after refinement)
      int NrNewNodes;
      /// Barycentric coordinates of the new nodes
      double* NewNodeCoords;
      /// The number of the new edges (after refinement)
      int NrNewEdges;
      /// Endpoints of the new edges
      int* NewEdges;
    };
    
  };
};

#endif
