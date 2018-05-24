#include "linalg/Vector.h"

namespace chemfe{
  namespace linalg{

    Vector::Vector(const size_t n) : n(n), data(std::vector<double>(n)) {}

    size_t Vector::size()
    {
      return n;
    }

    double& Vector::operator[](const size_t i)
    {
      return data.at(i);
    }
  }
}
