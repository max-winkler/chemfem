#include <algorithm>
#include <iterator>
#include <iomanip>

#include "linalg/SparseMatrix.h"

namespace chemfem{
  namespace linalg{

    SparseMatrix::SparseMatrix(size_t m, size_t n) : m(m), n(n), nnz(0) {}

    SparseMatrix::SparseMatrix(const SparseMatrix& M) : m(M.m), n(M.n), nnz(M.nnz)
    {
      Col = new size_t[nnz];
      Row = new size_t[m+1];
      Entry = new double[nnz];

      std::copy(M.Col, M.Col+nnz, Col);
      std::copy(M.Row, M.Row+m+1, Row);
      std::copy(M.Entry, M.Entry+nnz, Entry);
    }

    IdentityMatrix::IdentityMatrix(size_t m) : SparseMatrix(m, m)
    {
      nnz = m;
      Col = new size_t[m];
      Row = new size_t[m+1];
      Entry = new double[m];

      for(int i=0; i<m; ++i)
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
	  size_t cur_col;
	  for(size_t j=M.Row[i]; j<M.Row[i+1]; ++j)
	    {
	      cur_col = M.Col[j];
	      for(size_t k=last_col; k<cur_col; ++k)
		os << std::setw(5) << 0 << " ";

	      last_col = cur_col + 1;
	      
	      os << std::setw(5) << M.Entry[j] << " ";
	    }
	  for(size_t k=cur_col+1; k<M.n; ++k)
	    os << std::setw(5) << 0 << " ";
	  os << "]\n";
	}
      
      return os;
    }    
  };
};
