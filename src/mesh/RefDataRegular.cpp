#include "mesh/RefDataRegular.h"

static const int DatNewCells[][3] = {{0, 3, 5}, {1, 4, 3}, {2, 5, 4}, {3, 4, 5}};
static const int DatNewCellEdges[][3] = {{0, 5, 8}, {2, 6, 1}, {4, 7, 3}, {6, 7, 8}};
static const double DatNewNodeCoords[][2] = {{0.5, 0.0}, {0.5, 0.5}, {0.0, 0.5}};
static const int DatNewEdges[][2] = {{0, 3}, {3, 1}, {1, 4}, {4, 2}, {2, 5}, {5, 0},
				  {3, 4}, {4, 5}, {5, 3}};
static const int DatNewEdgeVertices[] = {3, 4, 5};

namespace chemfem{
  namespace mesh{

    RefDataRegular::RefDataRegular()
      : RefData(4, 6, 9)
    {            
      NewCells = DatNewCells;
      NewCellEdges = DatNewCellEdges;
      NewNodeCoords = DatNewNodeCoords;
      NewEdges = DatNewEdges;
      NewEdgeVertices = DatNewEdgeVertices;      
    }
    
  }
}
