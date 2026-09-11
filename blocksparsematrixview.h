#ifndef BLOCKSPARSEMATRIXVIEW_H
#define BLOCKSPARSEMATRIXVIEW_H

#include <vector>
#include "matrix.h"
#include "utils.hpp"
#include "blockallocator.h"
#include "blockallocator.hpp"

template<typename Num>
class BlockSparseMatrixView final{
  using Mat = Matrix<Num,BlockAllocator<Num>>;
  public:
    explicit BlockSparseMatrixView() = default;//empty view
    //full view onto a BlockSparseMatrix
    explicit BlockSparseMatrixView(BlockSparseMatrix<Num>& in);
    //sub-view onto a BlockSparseMatrix
    explicit BlockSparseMatrixView(BlockSparseMatrix<Num>& in, size_t row_block_start, size_t row_block_size, size_t row_block_stride, size_t col_block_start, size_t col_block_size, size_t col_block_stride);

    size_t nrowblocks() const {return _nrowblocks;}
    size_t ncolblocks() const {return _ncolblocks;}
    Mat*& block_ptr      (const size_t row_block, const size_t col_block)       & {return _block_ptrs[ij(row_block,col_block,_nrowblocks)];}
    const Mat*& block_ptr(const size_t row_block, const size_t col_block) const & {return _block_ptrs[ij(row_block,col_block,_nrowblocks)];}
    Mat&       block(const size_t row_block, const size_t col_block)       & {return *(this->block_ptr(row_block,col_block));}
    const Mat& block(const size_t row_block, const size_t col_block) const & {return *(this->block_ptr(row_block,col_block));}

    const BlockSparseMatrix<Num>* matrix() const {return _matrix;}
          BlockSparseMatrix<Num>* matrix()       {return _matrix;}

  private:
    //add offset/stride information as data members here?
    size_t _nrowblocks = 0;
    size_t _ncolblocks = 0;
    std::vector<Mat*> _block_ptrs;
    BlockSparseMatrix<Num>* _matrix;//ptr to the underlying (owning) matrix
};

//full view onto a BlockSparseMatrix (delegate to full constructor)
template<typename Num>
BlockSparseMatrixView<Num>::BlockSparseMatrixView(BlockSparseMatrix<Num>& in)
  : BlockSparseMatrixView(in,0,in.nrowblocks(),1,0,in.ncolblocks(),1) {}

template<typename Num>
BlockSparseMatrixView<Num>::BlockSparseMatrixView(BlockSparseMatrix<Num>& in, 
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
      this->block_ptr(rb_sub,cb_sub) = &in.block(rb_glob,cb_glob);
    }
  }
}


#endif

