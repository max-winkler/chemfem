#include "mesh/Cell.h"

namespace chemfem{
  namespace mesh{

    Cell::Cell(Node& n1, Node& n2, Node& n3)
    {
      LocNode[0] = &n1;
      LocNode[1] = &n2;
      LocNode[2] = &n3;
    }

    double Cell::Determinant() const
    {
      double x[3], y[3];
      for(int k=0; k<3; ++k)
	{
	  x[k] = LocNode[k]->getX();
	  y[k] = LocNode[k]->getY();
	}
      
      return (x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]);
    }

    using chemfem::linalg::DenseMatrix;
    
    DenseMatrix Cell::Jacobian() const
    {
      double x[3], y[3];
      for(int k=0; k<3; ++k)
	{
	  x[k] = LocNode[k]->getX();
	  y[k] = LocNode[k]->getY();
	}
      
      DenseMatrix Jac(2,2);
      Jac.data[0] =  x[1] - x[0];
      Jac.data[1] =  x[2] - x[0];
      Jac.data[2] =  y[1] - y[0];
      Jac.data[3] =  y[2] - y[0];

      return Jac;
    }

    size_t Cell::Index() const
    {
      return index;
    }

    void Cell::SetIndex(size_t index)
    {
      this->index = index;
    }

    int Cell::EdgeIndex(const Edge& edge) const
    {
      int edge_ind;
      for(edge_ind=0; edge_ind<3; ++edge_ind)
	if(LocEdge[edge_ind] == &edge) break;

      if(edge_ind == 3) edge_ind = -1;

      return edge_ind;
    }

    const Node& Cell::GetNode(int i) const
    {
      return *(LocNode[i]);
    }

    void Cell::SetRefType(RefType ref_type)
    {
      this->ref_type = ref_type;
    }
    
  };
};
