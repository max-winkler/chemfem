#include "linalg/Vector.h"

namespace chemfe{
  namespace linalg{

    Vector::Vector(const size_t n) : n(n)
    {
      data = new double[n];
    }

    Vector::~Vector()
    {
      delete[] data;
    }
    
    size_t Vector::size()
    {
      return n;
    }

    double& Vector::operator[](const size_t i)
    {
      return data[i];
    }
  }
}
