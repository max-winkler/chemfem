#ifndef _SPARSE_MATRIX_INSERTER_
#define _SPARSE_MATRIX_INSERTER_

#include <vector>

#include "linalg/SparseMatrix.h"

namespace chemfem{
  namespace linalg{
    class SparseMatrixInserter
    {      
    public:
      /**
       * Constructor, initialized with the matrix that should be inserted
       */
      SparseMatrixInserter(SparseMatrix&);

      /**
       * Build up the matrix
       */
      void Build();

      /**
       * Insert an entry in the specified row and column. Double entries are added up
       */
      void Insert(size_t, size_t, double);

    private:
      /**
       * Data structure for a single matrix entry. Used in the SparseMatrixInserter and is 
       * basically used for sorting entries in the Entries vector.
       */
      struct MatrixEntry
      {
      MatrixEntry(size_t row, size_t col, double entry) : col(col), row(row), entry(entry) {}
	size_t col, row;
	double entry;

	friend bool operator<(const MatrixEntry& A, const MatrixEntry& B)
	{
	  if(A.row < B.row) return true;
	  if(A.row > B.row) return false;
	  if(A.col < B.col) return true;
	  else
	    return false;
	}

	friend bool operator==(const MatrixEntry&A, const MatrixEntry& B)
	{
	  return A.row == B.row && A.col == B.col;
	}
      };
      
      std::vector<MatrixEntry> Entries;

      SparseMatrix& M;
    };
  };
};

#endif
