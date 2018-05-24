#include "mesh/UnitSquareMesh.h"

namespace chemfem{
  namespace mesh{

    UnitSquareMesh::UnitSquareMesh(size_t n)
    {
      // Node list
      for(size_t i=0; i<n; ++i)	
	for(size_t j=0; j<n; ++j)
	  {
	    Node node((double)i/n, (double)j/n);
	    Nodes.push_back(node);
	  }
      
      // Triangle list
      for(size_t i=0; i<n-1; ++i)
	for(size_t j=0; j<n-1; ++j)
	  {
	    Cells.push_back(Cell(i*n+j, i*n+j+1, (i+1)*n+j+1));
	    Cells.push_back(Cell(i*n+j, (i+1)*n+j+1, (i+1)*n+j));
	  }
    }
    
  }
}
