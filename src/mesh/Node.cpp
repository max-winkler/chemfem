#include "mesh/Node.h"

namespace chemfem{
  namespace mesh{

    Node::Node(size_t index, double x, double y) : index(index), x(x), y(y) {}

    double Node::getX() const { return x; }
    double Node::getY() const { return y; }

    size_t Node::Index() const
    {
      return index;
    }
    
  }
}
