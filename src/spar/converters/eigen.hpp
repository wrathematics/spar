// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_CONVERTERS_EIGEN_H
#define SPAR_CONVERTERS_EIGEN_H
#pragma once


#include <cstdint>

#include <Eigen/SparseCore>

#include "../arraytools/src/arraytools.hpp"


namespace spar
{
  template <typename INDEX, typename SCALAR>
  class dvec;

  template <typename INDEX, typename SCALAR>
  class spvec;

  template <typename INDEX, typename SCALAR>
  class spmat;
  
  namespace conv
  {
    /**
      @brief Shorthand for supported Eigen sparse matrix type.
      
      @tparam INDEX should be some kind of fundamental indexing type, like `int`
      or `uint16_t`.
      @tparam SCALAR should be a fundamental numeric type like `int` or `float`.
     */
    template <typename INDEX, typename SCALAR>
    using eigen_t = Eigen::SparseMatrix<SCALAR, Eigen::StorageOptions::ColMajor, INDEX>;
    
    
    
    /**
      @brief Convert an `spmat` object into an `Eigen::SparseMatrix` object.
      
      @param[in] x The input `spmat` object.
      
      @return The return sparse matrix.
      
      @allocs The return object is roughly of size:
      `sizeof(INDEX)*(2 + nnz + (n+1)) + sizeof(SCALAR)*nnz`.
      
      @except If a memory allocation fails, a `bad_alloc` exception will be
      thrown.
      
      @tparam INDEX should be some kind of fundamental indexing type, like `int`
      or `uint16_t`.
      @tparam SCALAR should be a fundamental numeric type like `int` or `float`.
     */
    template <typename INDEX, typename SCALAR>
    static inline eigen_t<INDEX, SCALAR> spmat_to_eigen(const spmat<INDEX, SCALAR> &x)
    {
      const INDEX m = x.nrows();
      const INDEX n = x.ncols();
      const INDEX nnz = x.get_nnz();
      
      eigen_t<INDEX, SCALAR> s(m, n);
      s.makeCompressed();
      s.resizeNonZeros(nnz);
      
      arraytools::copy(nnz, x.data_ptr(), s.valuePtr());
      arraytools::copy(nnz, x.index_ptr(), s.innerIndexPtr());
      arraytools::copy(n+1, x.col_ptr(), s.outerIndexPtr());
      
      return s;
    }
    
    
    
    /**
      @brief Convert an `Eigen::SparseMatrix` object into an `spmat` object.
      
      @param[in] s The input `Eigen::SparseMatrix` object.
      
      @return The return sparse matrix.
      
      @allocs The return object is roughly of size:
      `sizeof(INDEX)*(2 + nnz + (n+1)) + sizeof(SCALAR)*nnz`.
      
      @except If a memory allocation fails, a `bad_alloc` exception will be
      thrown.
      
      @tparam INDEX should be some kind of fundamental indexing type, like `int`
      or `uint16_t`.
      @tparam SCALAR should be a fundamental numeric type like `int` or `float`.
     */
    template <typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> eigen_to_spmat(const eigen_t<INDEX, SCALAR> &s)
    {
      const INDEX m = s.rows();
      const INDEX n = s.cols();
      const INDEX nnz = s.nonZeros();
      const INDEX p_len = s.outerSize();
      
      spmat<INDEX, SCALAR> x(m, n, nnz);
      
      arraytools::copy(nnz, s.valuePtr(), x.data_ptr());
      arraytools::copy(nnz, s.innerIndexPtr(), x.index_ptr());
      arraytools::copy(p_len, s.outerIndexPtr(), x.col_ptr());
      
      return x;
    }
    
    
    
    /**
      @brief Convert a column of an `Eigen::SparseMatrix` object into a dense
      `dvec` object.
      
      @param[in] col The desired column.
      @param[in] s The input `Eigen::SparseMatrix` object.
      
      @return The return dense vector.
      
      @allocs The return object is allocated to as many elements as the number
      of rows in the sparse matrix.
      
      @except If a memory allocation fails, a `bad_alloc` exception will be
      thrown.
      
      @tparam INDEX should be some kind of fundamental indexing type, like `int`
      or `uint16_t`.
      @tparam SCALAR should be a fundamental numeric type like `int` or `float`.
     */
    template <typename INDEX, typename SCALAR>
    static inline dvec<INDEX, SCALAR> eigen_to_dvec(const INDEX col, const eigen_t<INDEX, SCALAR> &s)
    {
      const INDEX m = s.rows();
      const INDEX n = s.cols();
      const INDEX nnz = s.nonZeros();
      
      const SCALAR *X = spmat_X(s);
      const INDEX *I = spmat_I(s);
      const INDEX *P = spmat_P(s);
      
      dvec<INDEX, SCALAR> d(m);
      
      const INDEX start = P[col];
      const INDEX end = (col == n) ? nnz : P[col+1];
      for (INDEX i=start; i<end; i++)
        d[ I[i] ] = X[i];
      
      return d;
    }
  }
  
  
  
  namespace internal
  {
    namespace get
    {
      template <typename INDEX, typename SCALAR>
      static inline void dim(const conv::eigen_t<INDEX, SCALAR> &x, INDEX *m, INDEX *n)
      {
        *m = (INDEX) x.rows();
        *n = (INDEX) x.cols();
      }
      
      
      
      template <typename INDEX, typename SCALAR>
      static inline void col(const INDEX j, const conv::eigen_t<INDEX, SCALAR> &x, spvec<INDEX, SCALAR> &s)
      {
        const INDEX *I = x.innerIndexPtr();
        const INDEX *P = x.outerIndexPtr();
        const SCALAR *X = x.valuePtr();
        
        const INDEX ind = P[j];
        if (P[j + 1] == ind)
        {
          s.zero();
          return;
        }
        
        const INDEX col_nnz = P[j + 1] - ind;
        s.set(col_nnz, I + ind, X + ind);
      }
      
      
      
      template <typename INDEX, typename SCALAR>
      static inline INDEX max_col_nnz(const conv::eigen_t<INDEX, SCALAR> &x)
      {
        const INDEX n = x.cols();
        const INDEX *P = x.outerIndexPtr();
        
        INDEX max_nnz = 0;
        for (INDEX col=0; col<n; col++)
        {
          INDEX col_nnz = P[col + 1] - P[col];
          if (col_nnz > max_nnz)
            max_nnz = col_nnz;
        }
        
        return max_nnz;
      }
    }
  }
}


#endif
