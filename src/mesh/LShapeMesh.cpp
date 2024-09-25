#include <algorithm>

#include "mesh/LShapeMesh.h"

namespace chemfem{
  namespace mesh{

    LShapeMesh::LShapeMesh(size_t n)
    {

      // Node list
      Nodes.reserve((2*n-1)*(2*n-1) - (n-1)*(n-1));

      size_t ctr = 0;
      
      for(size_t i=0; i<n-1; ++i)
        for(size_t j=0; j<n; ++j)
          Nodes.push_back(Node(ctr++, -1.+(double)j/(n-1), -1.+(double)i/(n-1)));

      for(size_t i=0; i<n; ++i)
        for(size_t j=0; j<2*n-1; ++j)
          Nodes.push_back(Node(ctr++, -1.+(double)j/(n-1), (double)i/(n-1)));
      
      // Triangle list
      Cells.reserve(6*(n-1)*(n-1));

      for(size_t i=0; i<n-1; ++i)
        for(size_t j=0; j<n-1; ++j)
          {
            Cells.push_back(Cell(i*n+j, i*n+j+1, (i+1)*n+j+1));
            Cells.push_back(Cell(i*n+j, (i+1)*n+j+1, (i+1)*n+j));
          }

      size_t s = n*(n-1);
      size_t r = 2*n-1;
      
      for(size_t i=0; i<n-1; ++i)
        for(size_t j=0; j<2*(n-1); ++j)
          {
            Cells.push_back(Cell(s+i*r+j, s+i*r+j+1, s+(i+1)*r+j+1));
            Cells.push_back(Cell(s+i*r+j, s+(i+1)*r+j+1, s+(i+1)*r+j));
          }
      
      CreateEdgeList();
    }    

  }
}
