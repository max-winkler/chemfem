#include "linalg/DenseMatrix.h"

namespace chemfem{
  namespace linalg{

    DenseMatrix::DenseMatrix(size_t m, size_t n) : m(m), n(n)
    {
      data = new double[m*n];
    }

    DenseMatrix::DenseMatrix(const DenseMatrix& Matrix)
    {
      copy(Matrix);
    }

    DenseMatrix& DenseMatrix::operator=(DenseMatrix& Matrix)
    {
      copy(Matrix);
      return *this;
    }

    void DenseMatrix::copy(const DenseMatrix& Matrix)
    {
      m = Matrix.m; n = Matrix.n;
      data = new double[m*n];
      std::copy(Matrix.data, Matrix.data+m*n, data);
    }

    DenseMatrix::~DenseMatrix()
    {
      delete[] data;
    }

    double DenseMatrix::Determinant() const
    {
      if(m != 2 || n != 2)
	{
	  std::cerr << "Implemented determinant routine only for 2x2 matrices.\n";
	  return 0;
	}
      
      return data[0]*data[3] - data[1]*data[2]; 
    }
    
    DenseMatrix DenseMatrix::Invert() const
    {
      DenseMatrix Inverse(m, m);
      
      if(m != 2 || n != 2)
	{
	  std::cerr << "Implemented inversion routine only for 2x2 matrices.\n";
	  return Inverse;
	}
      
      double Det = Determinant();
      
      Inverse.data[0] = data[3]/Det;
      Inverse.data[1] = -data[1]/Det;
      Inverse.data[2] = -data[2]/Det;
      Inverse.data[3] = data[0]/Det;

      return Inverse;
    }

    Vector DenseMatrix::operator*(const Vector& x) const
    {
      if(x.size() != n)
	{
	  std::cerr << "The matrix and the vector are incompatible for multiplication.\n";
	  return Vector(0);
	}

      Vector y(m);
      for(int i=0; i<m; ++i)
	{
	  y[i] = 0.;
	  for(int j=0; j<n; ++i)
	    y[i] += data[n*i+j]*x[j];
	}

      return y;
    }
  }
}
