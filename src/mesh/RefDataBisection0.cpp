#include "mesh/RefDataBisection0.h"

static const int DatNewCells[][3] = {{0, 1, 3}, {0, 3, 2}};
static const int DatNewCellEdges[][3] = {{0, 1, 4}, {4, 2, 3}};
static const double DatNewNodeCoords[][2] = {{0.5, 0.5}};
static const int DatNewEdges[][2] = {{0, 1}, {1, 3}, {3, 2}, {2, 0}, {3,0}};
static const int DatNewEdgeVertices[] = {-1, 3, -1};

namespace chemfem{
  namespace mesh{

    RefDataBisection0::RefDataBisection0()
      : RefData(2, 4, 5)
    {            
      NewCells = DatNewCells;
      NewCellEdges = DatNewCellEdges;
      NewNodeCoords = DatNewNodeCoords;
      NewEdges = DatNewEdges;
      NewEdgeVertices = DatNewEdgeVertices;      
    }
    
  }
}
