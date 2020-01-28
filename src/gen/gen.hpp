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
    
    
    template <typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> gen(const uint32_t seed, const INDEX nrows, const INDEX ncols)
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
    
    template <typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> gen(const INDEX nrows, const INDEX ncols)
    {
      const uint32_t seed = spar::internals::rand::get_seed();
      return gen<INDEX, SCALAR>(seed, nrows, ncols);
    }
  }
}


#endif
