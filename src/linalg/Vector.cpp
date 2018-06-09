#include <iostream>
#include <iomanip>
#include <math.h>

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
    const double& Vector::const_iterator::operator*() const
    {
      return *cur;
    }
        
    Vector::Vector() : n(0), data(NULL) {};
    
    Vector::Vector(const size_t n) : n(n)
    {
      data = new double[n]();      
    }

    Vector::Vector(const Vector& other) : n(other.n)
    {
      data = new double[n];
      std::copy(other.data, other.data+n, data);
    }
    
    Vector::~Vector()
    {
      if(data != NULL)
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
      if(data != NULL) // TODO: Reuse memory if possible
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

    Vector Vector::operator-(const Vector& b) const
    {
      Vector x(n);

      Vector::const_iterator it_a, it_b;
      size_t i;
      
      for(i=0, it_a = begin(), it_b = b.begin(); it_a != end(); ++it_a, ++it_b, ++i)
	x[i] = (*it_a) - (*it_b);

      return x;
    }
        
    void Vector::axpy(double a, const Vector& b, Vector& x) const
    {
      if(n != b.size())
	{
	  std::cerr << "Vectors are not compatible for summation.\n";
	  return;
	}

      Vector::const_iterator it_a, it_b;
      size_t i;
      
      for(i=0, it_a = begin(), it_b = b.begin(); it_a != end(); ++it_a, ++it_b, ++i)
	x[i] = a*(*it_a) + (*it_b);
    }

    Vector operator*(const double& a, const Vector& v)
    {
      Vector x(v.size());
      
      Vector::const_ierator it;
      int index;
      
      for(it = v.begin(), index=0; it != v.end(); ++it, ++index)
	x[index] = a*(*it);

      return x;
    }
    
    double Vector::Norm(double b) const
    {
      double norm = 0.;
      for(Vector::const_iterator it = begin(); it != end(); ++it)
	{
	  norm += pow(*it, b);
	}
      norm = pow(norm, 1./b);
    }
    
    Vector::const_iterator Vector::begin() const
    {
      return const_iterator(data);
    }

    Vector::const_iterator Vector::end() const
    {
      return const_iterator(data+n);
    }

    double dot(const Vector& a, const Vector& b)
    {
      double val = 0.;
      Vector::const_iterator it_a = a.begin();
      Vector::const_iterator it_b = b.begin();

      for(; it_a != a.end(), it_b != b.end(); ++it_a, ++it_b)
	val += (*it_a) * (*it_b);

      if(it_a != a.end() || it_b != b.end())
	{
	  std::cerr << "Error computing the dot product. Dimension of the vectors not equal.\n";
	  return 0;
	}
      return val;
    }

    std::ostream& operator<<(std::ostream& os, const Vector& Vec)
    {
      os << "[ ";
      for(Vector::const_iterator it = Vec.begin(); it != Vec.end(); ++it)	
	os << std::setw(5) << *it << " ";	  

      os << "]";
      return os;
    }
    
  }
}
