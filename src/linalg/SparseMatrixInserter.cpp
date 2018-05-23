#include <algorithm>

#include "linalg/SparseMatrixInserter.h"

namespace chemfem{
  namespace linalg{

    SparseMatrixInserter::SparseMatrixInserter(SparseMatrix& M) : M(M) {};

    void SparseMatrixInserter::insert(size_t i, size_t j, double a)
    {
      Entries.push_back(MatrixEntry(i, j, a));
    }
    
    void SparseMatrixInserter::build()
    {
      // Sort entries
      std::sort(Entries.begin(), Entries.end());

      // Determine number of non-zeros
      for(std::vector<MatrixEntry>::const_iterator it = Entries.begin();
	  it < Entries.end();)
	{
	  M.nnz++;
	  
	  std::vector<MatrixEntry>::const_iterator it_tmp = it;

	  // Sum up double entries
	  while(*it_tmp == *it) ++it;
	}

      M.Col = new size_t[M.nnz];
      size_t *RowFull = new size_t[M.nnz]; // Full row vector (uncompressed)
      M.Row = new size_t[M.m+1];
      M.Entry = new double[M.nnz];
      
      // Accummulate
      size_t entry_ctr = 0;
      for(std::vector<MatrixEntry>::const_iterator it = Entries.begin();
	  it < Entries.end();)
	{
	  double val = 0.;
	  std::vector<MatrixEntry>::const_iterator it_tmp = it;

	  // Sum up double entries
	  while(*it_tmp == *it)
	    {
	      val += it->entry;
	      ++it;
	    }

	  // Insert into Matrix
	  RowFull[entry_ctr] = it_tmp->row;
	  M.Col[entry_ctr] = it_tmp->col;
	  M.Entry[entry_ctr] = val;

	  ++entry_ctr;	 
	}

      // Compress row vector
      M.Row[0] = 0;

      size_t pos = 0, pos_tmp = 0;
      while(pos < M.nnz)
	{
	  size_t NrEntriesRow = 0;
	  while(RowFull[pos] == RowFull[pos+NrEntriesRow])
	    ++NrEntriesRow;
	  
	  M.Row[RowFull[pos]+1] = M.Row[RowFull[pos]] + NrEntriesRow;
	  pos = pos + NrEntriesRow;
	}
      delete[] RowFull;
    }    
  }
}
