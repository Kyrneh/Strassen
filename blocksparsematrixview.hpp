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

template<typename T1, typename T2>
typename T1::value_type dot(const BlockSparseMatrixView_basic<T1>& lhs,
                             const BlockSparseMatrixView_basic<T2>& rhs, const typename T1::value_type thresh){
  using Num = typename T1::value_type;//check that we cannot do float-double dots
  static_assert(std::is_same_v<Num,typename T2::value_type>,"dot: lhs and rhs must have the same Num type");
  //check that dimensions match
  assert(lhs.nblocks()    == rhs.nblocks());
  assert(lhs.nrowblocks() == rhs.nrowblocks());
  assert(lhs.ncolblocks() == rhs.ncolblocks());
  const size_t nijb = lhs.nblocks();
  const size_t njb  = lhs.ncolblocks();
  //adjust threshold to be per block instead of per element
  const Num thresh_per_block = thresh*static_cast<Num>(lhs.max_blocksize_row()*lhs.max_blocksize_col());

  Num retval = Num(0);
  //parallel loop over blocks
  #pragma omp parallel for schedule(dynamic) reduction(+: retval)
  for(size_t ijb=0;ijb<nijb;++ijb){
    const size_t ib = ijb/njb;
    const size_t jb = ijb%njb;
    const auto& lhs_block = lhs.block(ib,jb);
    const auto& rhs_block = rhs.block(ib,jb);
    if(lhs_block.size() != 0 && rhs_block.size() != 0){
      if(lhs_block.frobenius_norm()*rhs_block.frobenius_norm() >= thresh_per_block){
        retval += dot(lhs_block,rhs_block);
      }
    }
  }
  return retval;
}

#endif
