#include "mesh/Node.h"

namespace chemfem{
  namespace mesh{

    Node::Node(size_t Index, double x, double y) : Index(Index), x(x), y(y) {}

    double Node::getX() const { return x; }
    double Node::getY() const { return y; }
    
    
  }
}
