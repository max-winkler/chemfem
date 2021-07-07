#include "mesh/RefDataBisection2.h"

static const int DatNewCells[][3] = {{0, 3, 2}, {1, 2, 3}};
static const int DatNewCellEdges[][3] = {{0, 4, 3}, {2, 4, 1}};
static const double DatNewNodeCoords[][2] = {{0.5, 0.0}};
static const int DatNewEdges[][2] = {{0, 3}, {3, 1}, {1, 2}, {2, 0}, {2,3}};
static const int DatNewEdgeVertices[] = {3, -1, -1};

namespace chemfem{
  namespace mesh{

    RefDataBisection2::RefDataBisection2()
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
