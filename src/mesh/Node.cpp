#include "mesh/Node.h"

namespace chemfem{
  namespace mesh{

    Node::Node(double x, double y) : x(x), y(y) {}

    double Node::getX() const { return x; }
    double Node::getY() const { return y; }
    
    
  }
}
