// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_CONVERTERS_EIGEN_H
#define SPAR_CONVERTERS_EIGEN_H
#pragma once


#include <cstdint>

#include <Eigen/SparseCore>

#include "../arraytools/src/arraytools.hpp"
#include "../spar.hpp"


namespace spar
{
  namespace conv
  {
    namespace
    {
      template <typename INDEX, typename SCALAR>
      using eigen_t = Eigen::SparseMatrix<SCALAR, Eigen::StorageOptions::ColMajor, INDEX>;
    }
    
    
    
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
    
    
    
    template <typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> spmat_to_s4(const eigen_t<INDEX, SCALAR> &s)
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
}


#endif
