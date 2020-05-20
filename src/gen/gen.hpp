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
#include <vector>

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
      // sample a vector of length `len` with values coming from the interval
      // [0, n-1]
      template <typename INDEX>
      static std::vector<INDEX> res_sampler(const uint32_t seed, const INDEX n,
        const INDEX len, const bool zero_based_indices=true)
      {
        const INDEX offset = (INDEX) 1 - zero_based_indices;
        
        std::vector<INDEX> x(len);
        
        #pragma omp for simd
        for (INDEX i=0; i<len; i++)
          x[i] = i + offset;
        
        std::mt19937 mt(seed);
        
        for (INDEX i=len; i<n; i++)
        {
          std::uniform_int_distribution<INDEX> dist(0, i);
          INDEX j = dist(mt);
          
          if (j < len)
            x[j] = i + offset;
        }
        
        std::sort(x.begin(), x.end());
        
        return x;
      }
    }
    
    /**
      @brief Generate a sparse matrix whose entries are randomly
      TODO
      
      @param[in] seed Random seed.
      @param[in] prop_dense Proportion of density. Should be a number between
      0 and 1. Other numbers will be truncated to 0 (those less than 0) or 1
      (those greater than 1).
      @param[in] nrows,ncols The dimensions of the return.
      
      @return A sparse matrix.
      
      @allocs The return matrix is allocated/resized as necessary.
      
      @except If a memory allocation fails, a `bad_alloc` exception will be
      thrown.
     */
    template <typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> rand(const uint32_t seed,
      const float prop_dense, const INDEX nrows, const INDEX ncols)
    {
      const INDEX nnz = (INDEX) (std::min(std::max(prop_dense, 0.f), 1.f) * nrows*ncols);
      spmat<INDEX, SCALAR> x(nrows, ncols, nnz);
      if (nnz == 0)
        return x;
      
      INDEX *I = x.index_ptr();
      INDEX *P = x.col_ptr();
      SCALAR *X = x.data_ptr();
      
      const INDEX plen = ncols + 1;
      
      std::vector<INDEX> indices = internals::res_sampler(seed, nrows*ncols, nnz);
      INDEX col = 0;
      for (INDEX ind=0; ind<nnz; ind++)
      {
        INDEX i = indices[ind] % nrows;
        INDEX j = (indices[ind] - i) / nrows;
        
        I[ind] = i;
        X[ind] = 1;
        
        if (col == j)
        {
          P[col] = ind;
          col++;
        }
        
        while (col < j)
        {
          P[col] = ind;
          col++;
        }
      }
      
      for (INDEX ind=col; ind<plen; ind++)
        P[ind] += nnz;
      
      x.update_nnz();
      
      return x;
    }
    
    /// \overload
    template <typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> rand(const float prop_dense,
      const INDEX nrows, const INDEX ncols)
    {
      const uint32_t seed = spar::internals::rand::get_seed();
      return rand<INDEX, SCALAR>(seed, prop_dense, nrows, ncols);
    }
    
    
    
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
