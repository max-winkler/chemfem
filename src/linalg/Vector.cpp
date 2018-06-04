#include "linalg/Vector.h"

namespace chemfem{
  namespace linalg{

    Vector::const_iterator::const_iterator() : cur(NULL) {}
    Vector::const_iterator::const_iterator(const const_iterator& it) : cur(it.cur) {}
    Vector::const_iterator::const_iterator(const double* cur) : cur(cur) {}
    Vector::const_iterator& Vector::const_iterator::operator=(const const_iterator& it)
    {
      cur = it.cur;
      return *this;
    }
    bool Vector::const_iterator::operator==(const const_iterator& other)
    {
      return cur == other.cur;
    }
    bool Vector::const_iterator::operator!=(const const_iterator& other)
    {
      return cur != other.cur;
    }
    Vector::const_iterator& Vector::const_iterator::operator++()
    {
      ++cur;
    }
    const double Vector::const_iterator::operator*() const
    {
      return *cur;
    }
    
    
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

    const double& Vector::operator[](const size_t i) const
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

    Vector Vector::operator+(const Vector& b) const
    {
      Vector x(n);
      axpy(1, b, x);
      return x;
    }
        
    void Vector::axpy(double a, const Vector& b, Vector& x) const
    {
      if(n != b.size())
	{
	  std::cerr << "Vectors are not compatible for summation.\n";
	  return;
	}

      for(int i=0; i<n; ++i)
	x[i] = a*data[i] + b[i];
    }

    Vector::const_iterator Vector::begin() const
    {
      return const_iterator(data);
    }

    Vector::const_iterator Vector::end() const
    {
      return const_iterator(data+n);
    }
  }
}
