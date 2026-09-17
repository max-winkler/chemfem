#include <algorithm>
#include <iterator>
#include <iomanip>

#include "linalg/SparseMatrix.h"

namespace chemfem{
  namespace linalg{

    SparseMatrix::SparseMatrix(size_t m, size_t n) : m(m), n(n), nnz(0) {}

    SparseMatrix::SparseMatrix(const SparseMatrix& M) : m(M.m), n(M.n), nnz(M.nnz)
    {
      // A matrix that has not been built yet holds null pointers, copying from those would
      // read from nowhere
      if(M.Col && nnz)
	{
	  Col = new size_t[nnz];
	  Entry = new double[nnz];

	  std::copy(M.Col, M.Col+nnz, Col);
	  std::copy(M.Entry, M.Entry+nnz, Entry);
	}

      if(M.Row)
	{
	  Row = new size_t[m+1];
	  std::copy(M.Row, M.Row+m+1, Row);
	}
    }

    SparseMatrix::~SparseMatrix()
    {
      delete[] Col;
      delete[] Row;
      delete[] Entry;
    }

    SparseMatrix& SparseMatrix::operator=(const SparseMatrix& M)
    {
      if(this == &M)
	return *this;

      delete[] Col;
      delete[] Row;
      delete[] Entry;

      Col = NULL;
      Row = NULL;
      Entry = NULL;

      m = M.m;
      n = M.n;
      nnz = M.nnz;

      if(M.Col && nnz)
	{
	  Col = new size_t[nnz];
	  Entry = new double[nnz];

	  std::copy(M.Col, M.Col+nnz, Col);
	  std::copy(M.Entry, M.Entry+nnz, Entry);
	}

      if(M.Row)
	{
	  Row = new size_t[m+1];
	  std::copy(M.Row, M.Row+m+1, Row);
	}

      return *this;
    }

    SparseMatrix::SparseMatrix(SparseMatrix&& M) noexcept
      : Col(M.Col), Row(M.Row), Entry(M.Entry), m(M.m), n(M.n), nnz(M.nnz)
    {
      M.Col = NULL;
      M.Row = NULL;
      M.Entry = NULL;
      M.nnz = 0;
    }

    SparseMatrix& SparseMatrix::operator=(SparseMatrix&& M) noexcept
    {
      if(this == &M)
	return *this;

      delete[] Col;
      delete[] Row;
      delete[] Entry;

      Col = M.Col;
      Row = M.Row;
      Entry = M.Entry;

      m = M.m;
      n = M.n;
      nnz = M.nnz;

      M.Col = NULL;
      M.Row = NULL;
      M.Entry = NULL;
      M.nnz = 0;

      return *this;
    }

    Vector SparseMatrix::operator*(const Vector& x) const
    {
      Vector y(m);
      size_t *row;
      size_t* col;
      double* entry = Entry;
      size_t i=0;
      
      for(row = Row, col = Col; row!=Row+m; ++row, ++i)
	for(; col != Col+*(row+1); ++col, ++entry)	
	  y[i] += (*entry)*x[*col];

      return y;
    }
    
    void SparseMatrix::EliminateRowAndColumn(size_t k)
    {
      if(k >= m || k >= n)
	{
	  std::cerr << "Error: The index " << k << " is outside the matrix.\n";
	  return;
	}

      bool diagonal = false;

      for(size_t i=0; i<m; ++i)
	for(size_t j=Row[i]; j<Row[i+1]; ++j)
	  {
	    if(i == k && Col[j] == k)
	      {
		Entry[j] = 1.;
		diagonal = true;
	      }
	    else if(i == k || Col[j] == k)
	      Entry[j] = 0.;
	  }

      if(!diagonal)
	std::cerr << "Error: The diagonal entry " << k << " is not in the sparsity pattern, "
		  << "the matrix is singular now.\n";
    }

    IdentityMatrix::IdentityMatrix(size_t m) : SparseMatrix(m, m)
    {
      nnz = m;
      Col = new size_t[m];
      Row = new size_t[m+1];
      Entry = new double[m];

      for(size_t i=0; i<m; ++i)
	{
	  Col[i] = i;
	  Row[i] = i;
	  Entry[i] = 1.;
	}
      Row[m] = m;
    }    
    
    std::ostream& operator<<(std::ostream& os, const SparseMatrix& M)
    {
      for(size_t i=0; i<M.m; ++i)
	{
	  os << " [ ";
	  
	  size_t last_col = 0;
	  size_t cur_col = 0;
	  for(size_t j=M.Row[i]; j<M.Row[i+1]; ++j)
	    {
	      cur_col = M.Col[j];
	      for(size_t k=last_col; k<cur_col; ++k)
		os << std::setw(10) << 0 << " ";

	      last_col = cur_col + 1;
	      
	      os << std::setw(10) << M.Entry[j] << " ";
	    }
	  for(size_t k=cur_col+1; k<M.n; ++k)
	    os << std::setw(10) << 0 << " ";
	  os << "]\n";
	}
      
      return os;
    }    
  };
};
