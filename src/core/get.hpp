// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_CORE_GET_H
#define SPAR_CORE_GET_H
#pragma once


template <typename INDEX, typename SCALAR>
class spvec;

template <typename INDEX, typename SCALAR>
class spmat;

namespace spar
{
  namespace internal
  {
    namespace get
    {
      template <typename INDEX, typename SCALAR>
      static inline void dim(const spmat<INDEX, SCALAR> &x, INDEX *m, INDEX *n)
      {
        *m = x.nrows();
        *n = x.ncols();
      }
      
      
      
      template <typename INDEX, typename SCALAR>
      static inline void col(const INDEX j, const spmat<INDEX, SCALAR> &x, spvec<INDEX, SCALAR> &s)
      {
        x.get_col(j, s);
      }
      
      
      
      template <typename INDEX, typename SCALAR>
      static inline INDEX max_col_nnz(const spmat<INDEX, SCALAR> &x)
      {
        INDEX max_nnz = 0;
        
        const INDEX *P = x.col_ptr();
        const INDEX n = x.ncols();
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
