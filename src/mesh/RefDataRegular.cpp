#include "mesh/RefDataRegular.h"

namespace chemfem{
  namespace mesh{

    static const int NewCells_[][3] = {{0, 3, 5}, {1, 4, 3}, {2, 5, 4}, {3, 4, 5}};
    static const double NewNodeCoords_[][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1},
					       {0.5, 0.5, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5}};
    static const int NewEdges_[][2] = {{0,3}, {3,1}, {1,4}, {4,2}, {2,5}, {5,0},
				       {3,4}, {4,5}, {5,3}};
    
    RefDataRegular::RefDataRegular()
      : RefData(NrNewCells, NrNewNodes, NrNewEdges)
    {
      NewCells = (const int*) NewCells_;
      NewNodeCoords = (const double*) NewNodeCoords_;
      NewEdges = (const int*) NewEdges_;
    }
    
  }
}
