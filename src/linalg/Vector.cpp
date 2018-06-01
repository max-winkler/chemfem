#include "linalg/Vector.h"

namespace chemfem{
  namespace linalg{

    Vector::Vector() : n(0), data(NULL) {};
    
    Vector::Vector(const size_t n) : n(n)
    {
      data = new double[n];
    }

    Vector::~Vector()
    {
      delete[] data;
    }
    
    size_t Vector::size() const
    {
      return n;
    }

    double& Vector::operator[](const size_t i)
    {
      return data[i];
    }

    Vector& Vector::operator=(const Vector& v)
    {
      size_t n_ = v.size();

      n = n_;
      if(data != NULL)
	delete[] data;
      data = new double[n];
      
      std::copy(v.data, v.data+n, data);
      return *this;
    }

    Vector& Vector::operator=(const double* data_in)
    {
      std::copy(data_in, data_in+n, data);
    }
  }
}
