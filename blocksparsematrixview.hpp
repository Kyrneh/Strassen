#ifndef BLOCKSPARSEMATRIXVIEW_HPP
#define BLOCKSPARSEMATRIXVIEW_HPP
#include "blocksparsematrixview.h"
#include "blocksparsematrix.h"

//full view onto a BlockSparseMatrix (delegate to full constructor)
template<typename T>
BlockSparseMatrixView_basic<T>::BlockSparseMatrixView_basic(T& in)
  : BlockSparseMatrixView_basic(in,0,in.nrowblocks(),1,0,in.ncolblocks(),1) {}

template<typename T>
BlockSparseMatrixView_basic<T>::BlockSparseMatrixView_basic(T& in,
    const size_t row_block_start, const size_t row_block_size, const size_t row_block_stride,
    const size_t col_block_start, const size_t col_block_size, const size_t col_block_stride)
  : _nrowblocks(row_block_size),
    _ncolblocks(col_block_size),
    _block_ptrs(row_block_size*col_block_size,nullptr),
    _matrix(&in)
{
  //create a pointer to each element in the sub-block
  for(size_t rb_sub=0;rb_sub < row_block_size; ++rb_sub){
    const size_t rb_glob = row_block_start + rb_sub*row_block_stride;
    for(size_t cb_sub=0;cb_sub < col_block_size; ++cb_sub){
      const size_t cb_glob = col_block_start + cb_sub*col_block_stride;
      _block_ptrs[ij(rb_sub,cb_sub,_nrowblocks)] = &in.block(rb_glob,cb_glob);
    }
  }
}

template<typename T>
template<typename U, typename>
BlockSparseMatrixView_basic<T>::BlockSparseMatrixView_basic(const BlockSparseMatrixView_basic<U>& in)
  : _nrowblocks(in.nrowblocks()),
    _ncolblocks(in.ncolblocks()),
    _block_ptrs(in.nblocks()),
    _matrix(in.matrix())
{
  for(size_t rb=0;rb<_nrowblocks;++rb)
    for(size_t cb=0;cb<_ncolblocks;++cb)
      _block_ptrs[ij(rb,cb,_nrowblocks)] = in.block_ptr(rb,cb);//Mat* -> const Mat*, always valid
}

template<typename T>
size_t BlockSparseMatrixView_basic<T>::max_blocksize_row() const {return _matrix->max_blocksize_row();}
template<typename T>
size_t BlockSparseMatrixView_basic<T>::max_blocksize_col() const {return _matrix->max_blocksize_col();}

#endif
