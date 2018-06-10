#include "mesh/UnitSquareMesh.h"

namespace chemfem{
  namespace mesh{

    UnitSquareMesh::UnitSquareMesh(size_t n)
    {
      // Node list
      Nodes.reserve(n*n);				
      for(size_t i=0; i<n; ++i)
	for(size_t j=0; j<n; ++j)
	  Nodes.push_back(Node(i*n+j, (double)j/(n-1), (double)i/(n-1)));
      
      // Triangle list
      Cells.reserve(2*(n-1)*(n-1));
      for(size_t i=0; i<n-1; ++i)
	for(size_t j=0; j<n-1; ++j)
	  {
	    //Cells.push_back(Cell(i*n+j, i*n+j+1, (i+1)*n+j+1));
	    //Cells.push_back(Cell(i*n+j, (i+1)*n+j+1, (i+1)*n+j));
	    
	    Cells.push_back(Cell(Nodes[i*n+j], Nodes[i*n+j+1], Nodes[(i+1)*n+j]));
	    Cells.push_back(Cell(Nodes[i*n+j+1], Nodes[(i+1)*n+j+1], Nodes[(i+1)*n+j]));
	  }

      CreateEdgeList();
    }
    
  }
}
