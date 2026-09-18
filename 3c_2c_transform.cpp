#include "matrix.h"
#include "matrix.hpp"
#include "blocksparsematrix.h"
#include "blocksparsematrix.hpp"
#include "strassen.hpp"
#include <chrono>
#include <cstdio>

int main(int argc, char** argv){
  if(argc != 8){
    puts("Usage: <program> <N_aux> <N_vec2> <3c filename> <2c filename> <blocksize> <threshold> <chunk_size>");
    return 1;
  }
  const size_t N_aux = (size_t)std::stol(argv[1]);
  const size_t N_vec2 = (size_t)std::stol(argv[2]);
  const std::string fn_3c = std::string(argv[3]);
  const std::string fn_2c = std::string(argv[4]);
  const size_t bs = (size_t)std::stol(argv[5]);
  const double thresh_mult = std::stod(argv[6]);
  const size_t chunk_size = std::stol(argv[7]);
  if(chunk_size > N_aux){ 
    puts("chunk_size needs to be <= N_aux!");
    return(1);
  }
  if(chunk_size > N_vec2){ 
    puts("chunk_size needs to be <= N_vec2!");
    return(1);
  }

  Matrix<double> ints_3c(N_vec2, N_aux);
  ints_3c.read_from_file(fn_3c.c_str());
  Matrix<double> ints_2c(N_aux, N_aux);
  ints_2c.read_from_file(fn_2c.c_str());

  //(mu nu |P)^' = \sum_Q (mu nu|Q) (Q|P)^{-1/2}
  Matrix<double> ints_3c_transformed(N_vec2, N_aux);

  printf("--- dense matmult  ---\n");
  for(int i=0;i<1;++i)
  {
    const auto start=std::chrono::steady_clock::now();
    matmult(ints_3c_transformed,ints_3c,false,ints_2c,false,1.0,0.0);
    const auto end=std::chrono::steady_clock::now();
    const double us=(double)std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
    printf("  [%d] %.4f s  (%.4f GFLOPs)\n",i+1,1e-6*us,2e-3*(double)(N_vec2*N_aux*N_aux)/us);
  }
  const double L2_norm_of_output = ints_3c_transformed.calc_frobenius_norm();
  printf("L2 norm of output = %e\n",L2_norm_of_output);


#if 1
  {
    BlockSparseMatrix<double> ints_3c_bs(ints_3c,bs,bs,0.0);
    //ints_3c = Matrix<double>();
    BlockSparseMatrix<double> ints_2c_bs(ints_2c,bs,bs,0.0);
    //ints_2c = Matrix<double>();
    BlockSparseMatrix<double> ints_3c_transformed_bs(N_vec2, N_aux,bs,bs,0.0);

    printf("\n--- BSM matmult (thresh=%.2e) ---\n",thresh_mult);

    for(int i=0;i<1;++i)
    {
      ints_3c_transformed_bs.fill_with_values(0.e0);
      const auto start=std::chrono::steady_clock::now();
      matmult(ints_3c_transformed_bs,ints_3c_bs,false,ints_2c_bs,false,thresh_mult,1.0,1.0);
      const auto end=std::chrono::steady_clock::now();
      const double us=(double)std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
      printf("  [%d] %.4f s  (%.4f GFLOPs)",i+1,1e-6*us,2e-3*(double)(N_vec2*N_aux*N_aux)/us);
      printf(" relative RMSD = %e\n",(ints_3c_transformed_bs.to_matrix()-ints_3c_transformed).calc_frobenius_norm()/L2_norm_of_output);
    }
  }
#endif

#if 0
  {
    printf("\n--- Strassen matmult (thresh=%.2e) ---\n",thresh_mult);
    //cut the molecule into perfectly sized chunks and then perform Strassen mults on those chunks in parallel
    //const size_t chunk_size = bs*(size_t)(std::exp2(std::floor(std::log2((double)(N_aux/bs)+0.5)))+0.5);
    printf("chunk_size = %lu\n",chunk_size);
    Matrix<double>& ints_3c_transformed_strassen = ints_3c_transformed;

    Matrix<double> ints_2c_chopped(chunk_size,chunk_size);
    #pragma omp for parallel schedule(static) collapse(2)
    for (size_t col=0;col<chunk_size;++col){
      for (size_t row=0;row<chunk_size;++row){
        assert(row < ints_2c.nrow());
        assert(col < ints_2c.ncol());
        ints_2c_chopped.elem(row,col) = ints_2c.elem(row,col);
      }
    }
    const size_t n_chunks = N_vec2/chunk_size;

    for(int i=0;i<1;++i)
    {
      const auto start=std::chrono::steady_clock::now();
      #pragma omp parallel
      {
        BlockSparseMatrix<double> ints_2c_bs(ints_2c_chopped,bs,bs,0.0);
        Matrix<double> ints_3c_chopped(0.e0,chunk_size,chunk_size);
        //Matrix<double> ints_3c_transformed_chopped(0.e0,chunk_size,chunk_size);
        BlockSparseMatrix<double> ints_3c_bs(ints_3c_chopped,bs,bs,0.0);

        BlockSparseMatrix<double> ints_3c_transformed_bs(ints_3c_chopped,bs,bs,0.0);
        #pragma omp for schedule(static)
        for (size_t chunk=0;chunk<n_chunks; ++chunk){
          //create the chunk_size*chunk_size bsmat for 3c ints
          for (size_t col=0;col<chunk_size;++col){
            for (size_t row=0;row<chunk_size;++row){
              assert(row < ints_3c.nrow());
              assert(chunk*chunk_size+row < ints_3c.nrow());
              assert(col < ints_3c.ncol());
              ints_3c_chopped.elem(row,col) = ints_3c.elem(chunk*chunk_size+row,col);
            }
          }
          ints_3c_bs.copy_from_input_matrix(ints_3c_chopped);
          //ints_3c_transformed_bs.fill_with_values(0.e0);
          //do the the transformation (matmult)
          //matmult(ints_3c_transformed_chopped,ints_3c_chopped,false,ints_2c_chopped,false,1.0,0.0);
          matmult_strassen_sparse(ints_3c_transformed_bs,ints_3c_bs,false,ints_2c_bs,false,thresh_mult,1.0,1.0);
          //matmult                (ints_3c_transformed_bs,ints_3c_bs,false,ints_2c_bs,false,thresh_mult,1.0,1.0);
          //transform back to dense matrix (reuse ints_3c_chopped buffer)
          ints_3c_transformed_bs.to_pointer(ints_3c_chopped.data_ptr());
          //write back to output array
          for (size_t col=0;col<chunk_size;++col){
            for (size_t row=0;row<chunk_size;++row){
              assert(row < ints_3c_transformed.nrow());
              assert(chunk*chunk_size+row < ints_3c_transformed.nrow());
              assert(col < ints_3c_transformed.ncol());
              ints_3c_transformed_strassen.elem(chunk*chunk_size+row,col) = ints_3c_chopped.elem(row,col);
              //ints_3c_transformed_strassen.elem(chunk*chunk_size+row,col) = ints_3c_transformed_chopped.elem(row,col);
            }
          }
        }//end omp for
      }//end omp parallel
      const auto end=std::chrono::steady_clock::now();
      const double us=(double)std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
      printf("  [%d] %.4f s  (%.4f GFLOPs)\n",i+1,1e-6*us,2e-3*(double)(n_chunks*chunk_size*chunk_size*chunk_size)/us);
      //printf(" relative RMSD = %e\n",(ints_3c_transformed_strassen-ints_3c_transformed).calc_frobenius_norm()/L2_norm_of_output);
    }//end timing loop
  }//end Strassen block
#endif

}

