// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_GEN_GEN_H
#define SPAR_GEN_GEN_H
#pragma once


#include <algorithm>
#include <cstdint>
#include <random>
#include <stdexcept>

#include "rand.hpp"

template <typename INDEX, typename SCALAR>
class spvec;

template <typename INDEX, typename SCALAR>
class spmat;

namespace spar
{
  /// @brief Random generators.
  namespace gen
  {
    namespace internals
    {
      template <typename INDEX>
      static inline float fudge(const INDEX i, const INDEX j, const INDEX n)
      {
        float a = std::abs((float)(i-j)/n);
        return std::max(1.f - 2.f*a - ((float)2/n), 0.f);
      }
    }
    
    
    
    /**
      @brief Generate a banded-like matrix. Non-zero values are generated at
      random with decreasing probability the farther from the diagonal they are.
      
      @param[in] seed Random seed.
      @param[in] nrows,ncols The dimensions of the return.
      
      @return A sparse matrix.
      
      @allocs The return matrix is allocated/resized as necessary.
      
      @except If a memory allocation fails, a `bad_alloc` exception will be
      thrown.
     */
    template <typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> bandish(const uint32_t seed,
      const INDEX nrows, const INDEX ncols)
    {
      const INDEX slen = (INDEX) nrows/2;
      spvec<INDEX, SCALAR> s(slen);
      
      const INDEX xlen = (INDEX) slen * std::min(nrows, ncols);
      spmat<INDEX, SCALAR> x(nrows, ncols, xlen);
      
      std::mt19937 mt(seed);
      for (INDEX j=0; j<ncols; j++)
      {
        s.zero();
        
        for (INDEX i=0; i<nrows; i++)
        {
          const float p = internals::fudge(i, j, nrows);
          std::binomial_distribution<int> dist(1, p);
          
          SCALAR d = dist(mt);
          if (d > 0)
            s.insert(i, 1);
        }
        
        x.insert(j, s);
      }
      
      return x;
    }
    
    /// \overload
    template <typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> bandish(const INDEX nrows,
      const INDEX ncols)
    {
      const uint32_t seed = spar::internals::rand::get_seed();
      return bandish<INDEX, SCALAR>(seed, nrows, ncols);
    }
    
    
    
    namespace internals
    {
      template <typename INDEX>
      static inline INDEX ntri(const INDEX n)
      {
        return (n*(n+1))/2;
      }
      
      template <typename INDEX>
      static inline INDEX banded_nnz(const INDEX band, const INDEX nrows,
        const INDEX ncols)
      {
        const INDEX m = nrows>ncols?nrows:ncols;
        const INDEX n = nrows>ncols?ncols:nrows;
        
        INDEX nnz = n + ntri(n-1)-ntri(n-band) + ntri(m-1)-ntri(m-band) - (m-1-n);
        if (m == n)
          nnz -= (m-1 - n);
        
        return nnz;
      }
    }
    
    /**
      @brief Generate a banded matrix. Non-zero values are generated at
      random with decreasing probability the farther from the diagonal they are.
      
      @param[in] band Size of the band.
      @param[in] nrows,ncols The dimensions of the return.
      
      @return A sparse matrix.
      
      @allocs The return matrix is allocated/resized as necessary.
      
      @except If a memory allocation fails, a `bad_alloc` exception will be
      thrown.
     */
    template <typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> banded(const INDEX band,
      const INDEX nrows, const INDEX ncols)
    {
      const INDEX nnz = internals::banded_nnz(band, nrows, ncols);
      spmat<INDEX, SCALAR> x(nrows, ncols, nnz);
      
      INDEX *I = x.index_ptr();
      INDEX *P = x.col_ptr();
      SCALAR *X = x.data_ptr();
      
      INDEX pos = 0;
      INDEX col_pos = 1;
      INDEX i = 0, j = 0;
      
      P[0] = 0;
      
      while (i < nrows)
      {
        if (i == j || (i < j && j < i + band) || (i > j && i < j + band))
        {
          I[pos] = i;
          X[pos] = 1;
          pos++;
          i++;
        }
        else if (i >= band + j)
        {
          i = (i>band?i-band:0);
          j++;
          
          P[col_pos] = pos;
          col_pos++;
          
          if (j >= ncols)
            break;
        }
        else
          i++;
      }
      
      x.update_nnz();
      return x;
    }
  }
}


#endif
