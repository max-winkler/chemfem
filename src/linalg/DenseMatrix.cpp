#include <iostream>
#include <iomanip>

#include "linalg/DenseMatrix.h"

namespace chemfem{
  namespace linalg{

    DenseMatrix::DenseMatrix(size_t m, size_t n) : m(m), n(n)
    {
      data = new double[m*n];
      std::fill(data, data+m*n, 0.);
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

    void DenseMatrix::Init(const double* entries)
    {
      std::copy(entries, entries + m*n, data);
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

    DenseMatrix DenseMatrix::Transpose() const
    {
      DenseMatrix Transposed(n, m);

      for(size_t i=0; i<n; ++i)
	for(size_t j=0; j<m; ++j)
	  Transposed.data[i*n+j] = data[j*m+i];
      
      return Transposed;
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
	  for(int j=0; j<n; ++j)
	    y[i] += data[n*i+j]*x[j];
	}

      return y;
    }

    std::ostream& operator<<(std::ostream& os, const DenseMatrix& M)
    {
      for(size_t i=0; i<M.m; ++i)
	{
	  os << " [ ";
	  	
	  for(size_t j=0; j<M.n; ++j)
	    os << std::setw(5) << M[i][j] << " ";

	  os << "]\n";
	}
      
      return os;
    }    
    
    DenseMatrix::DenseMatrixRow DenseMatrix::operator[](const size_t i)
    {
      return DenseMatrixRow(&(data[n*i]));
    }

    const DenseMatrix::DenseMatrixRow DenseMatrix::operator[](const size_t i) const
    {
      return DenseMatrixRow(&(data[n*i]));
    }

    DenseMatrix::DenseMatrixRow::DenseMatrixRow(double* RowPtr) : RowPtr(RowPtr) {}

    double& DenseMatrix::DenseMatrixRow::operator[](const size_t j)
    {
      return RowPtr[j];
    }

    const double& DenseMatrix::DenseMatrixRow::operator[](const size_t j) const
    {
      return RowPtr[j];
    }
  }
}
