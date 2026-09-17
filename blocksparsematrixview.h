#ifndef BLOCKSPARSEMATRIXVIEW_H
#define BLOCKSPARSEMATRIXVIEW_H

#include <vector>
#include <type_traits>
#include "matrix.h"
#include "utils.hpp"
#include "blockallocator.h"
#include "blockallocator.hpp"

//forward declaration: the (owning) matrix this view refers to
template<typename Num>
class BlockSparseMatrix;

//non-owning, span-like view onto (a sub-block-grid of) a BlockSparseMatrix<Num>.
//T is either BlockSparseMatrix<Num> (mutable view) or const BlockSparseMatrix<Num>
//(read-only view) -- mirrors std::span<T> vs std::span<const T>.
template<typename T>
class BlockSparseMatrixView_basic final{
    using Num = typename T::value_type;
    using Mat = Matrix<Num,BlockAllocator<Num>>;
    //block type as seen through this view: const-qualified iff T is const
    using ViewedMat = std::conditional_t<std::is_const_v<T>,const Mat,Mat>;
  public:
    explicit BlockSparseMatrixView_basic() = default;//empty view
    //full view onto a BlockSparseMatrix
    explicit BlockSparseMatrixView_basic(T& in);
    //sub-view onto a BlockSparseMatrix
    explicit BlockSparseMatrixView_basic(T& in, size_t row_block_start, size_t row_block_size, size_t row_block_stride, size_t col_block_start, size_t col_block_size, size_t col_block_stride);
    //implicit const-adding conversion (mutable view -> read-only view), mirrors span<T> -> span<const T>
    template<typename U, typename = std::enable_if_t<std::is_same_v<T,std::add_const_t<U>>>>
    BlockSparseMatrixView_basic(const BlockSparseMatrixView_basic<U>& in);

    //accessors to properties of the underlying (owning) matrix
    size_t nrowblocks() const {return _nrowblocks;}
    size_t ncolblocks() const {return _ncolblocks;}
    size_t nblocks()    const {return _block_ptrs.size();}
    size_t max_blocksize_row() const;
    size_t max_blocksize_col() const;

    //shallow-const accessors: constness of the result depends only on T, never on
    //the constness of *this (same as std::span::data()/operator[])
    ViewedMat* block_ptr(const size_t row_block, const size_t col_block) const {return _block_ptrs[ij(row_block,col_block,_nrowblocks)];}
    ViewedMat& block(const size_t row_block, const size_t col_block)     const {return *(this->block_ptr(row_block,col_block));}
    T* matrix() const {return _matrix;}

  private:
    //add offset/stride information as data members here?
    size_t _nrowblocks = 0;
    size_t _ncolblocks = 0;
    std::vector<ViewedMat*> _block_ptrs;
    T* _matrix = nullptr;//ptr to the underlying (owning) matrix
};

template<typename Num>
using BlockSparseMatrixView      = BlockSparseMatrixView_basic<BlockSparseMatrix<Num>>;
template<typename Num>
using ConstBlockSparseMatrixView = BlockSparseMatrixView_basic<const BlockSparseMatrix<Num>>;

#endif
