#include <cmath>

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

    double Node::Dist(const Node& other) const
    {
      return sqrt(pow(x-other.x, 2.) + pow(y-other.y, 2.));
    }
    
    std::ostream& operator<< (std::ostream& stream, const Node& node)
    {
      stream << "(" << node.x << "," << node.y << ")";
      return stream;
    }
    
    
  }
}
