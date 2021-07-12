#ifndef _REF_DATA_H_
#define _REF_DATA_H_

#include <cstring>

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
      friend class Mesh;
    public:

      /// Constructor which initailized the number of new cells, nodes and edges
      RefData(int, int, int);
      ~RefData();
      
      /// Returns the number of nodes
      int GetNrNodes();
      /// Returns the number of new edges
      int GetNrEdges();
      /// Returns number of edges that will be refined
      int GetNrRefinedEdges();
      /// Returns the number of cells         
      int GetNrCells();
      /// Returns a pointer to the coordinates of the i-th new node
      const double* GetNodeCoords(int);
      /// Returns a pointer to the endpoints of the edges (dimension is 2)
      const int* GetEdge(int);
      /// Returns a pointer to the vertices of the new cells (dimension is 3)
      const int* GetCell(int);
      /// Returns a pointer to the indices of the edges of the new cells (dimension is 3)
      const int* GetCellEdges(int);      
      /// Returns the index of the new vertex on the refined edges
      const int GetEdgeVertex(int);

    protected:
      /// The number of new cells (after refinement)
      int NrNewCells;
      /// The indices of the vertices of new cells
      const int (*NewCells)[3];
      /// The indices of the edges of new cells
      const int (*NewCellEdges)[3];            
      /// The number of new nodes (after refinement)
      int NrNewNodes;
      /// Barycentric coordinates of the new nodes
      const double (*NewNodeCoords)[2];
      /// The number of the new edges (after refinement)
      int NrNewEdges;
      /// Endpoints of the new edges
      const int (*NewEdges)[2];
      /// New vertices at the midpoint of an edge (-1 if no new vertex is inserted)
      const int *NewEdgeVertices;
    };
    
  };
};

#endif
