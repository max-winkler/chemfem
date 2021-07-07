#include "mesh/RefDataBisection1.h"

static const int DatNewCells[][3] = {{0, 1, 3}, {1, 2, 3}};
static const int DatNewCellEdges[][3] = {{0, 4, 3}, {1, 2, 4}};
static const double DatNewNodeCoords[][2] = {{0.0, 0.5}};
static const int DatNewEdges[][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {1,3}};
static const int DatNewEdgeVertices[] = {-1, -1, 3};

namespace chemfem{
  namespace mesh{

    RefDataBisection1::RefDataBisection1()
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
