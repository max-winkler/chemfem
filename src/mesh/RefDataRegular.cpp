#include "mesh/RefDataRegular.h"

namespace chemfem{
  namespace mesh{

    RefDataRegular::RefDataRegular()
      : RefData(4, 6, 3)
    {
      const int NewCells[][3] = {{0, 3, 5}, {1, 4, 3}, {2, 5, 4}, {3, 4, 5}};
      const double NewNodeCoords[][2] = {{0.5, 0.0}, {0.5, 0.5}, {0.0, 0.5}};
      const int NewEdges[][2] = {{0, 3}, {3, 1}, {1, 4}, {4, 2}, {2, 5}, {5, 0},
				 {3, 4}, {4, 5}, {5, 6}};

      this->NewCells = NewCells;
      this->NewNodeCoords = NewNodeCoords;
      this->NewEdges = NewEdges;
    }
    
  }
}
