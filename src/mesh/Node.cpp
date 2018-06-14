#include "mesh/Node.h"

namespace chemfem{
  namespace mesh{

    Node::Node(size_t index, double x, double y) : x(x),  y(y), index(index) {}

    double Node::getX() const { return x; }
    double Node::getY() const { return y; }

    size_t Node::Index() const
    {
      return index;
    }
    
  }
}
