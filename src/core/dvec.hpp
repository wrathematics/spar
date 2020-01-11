// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_CORE_DVEC_H
#define SPAR_CORE_DVEC_H
#pragma once


#include <iostream>

#include "../arraytools/src/arraytools.hpp"


template <typename INDEX, typename SCALAR>
class spvec;

template <typename INDEX, typename SCALAR>
class dvec
{
  public:
    dvec();
    dvec(INDEX len_);
    ~dvec();
    
    void resize(INDEX len_);
    void zero();
    SCALAR sum() const;
    void insert(const INDEX i, const SCALAR s);
    void update_nnz();
    template <typename INDEX_SRC, typename SCALAR_SRC>
    void set(const INDEX_SRC nnz_, const INDEX_SRC *I_, const SCALAR_SRC *X_);
    void set(const spvec<INDEX, SCALAR> &x);
    
    void print() const;
    
    const SCALAR operator[](INDEX i) const; // getter
    SCALAR& operator[](INDEX i); // setter
    
    INDEX get_nnz() const {return nnz;};
    INDEX get_len() const {return len;};
    SCALAR* data_ptr() {return X;};
    SCALAR* data_ptr() const {return X;};
  
  protected:
    INDEX nnz;
    INDEX len;
    SCALAR *X;
  
  private:
    void cleanup();
};



// ----------------------------------------------------------------------------
// constructor/destructor
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
dvec<INDEX, SCALAR>::dvec()
{
  X = NULL;
  
  nnz = 0;
  len = 0;
}



template <typename INDEX, typename SCALAR>
dvec<INDEX, SCALAR>::dvec(INDEX len_)
{
  arraytools::zero_alloc(len_, &X);
  arraytools::check_allocs(X);
  
  nnz = 0;
  len = len_;
}



template <typename INDEX, typename SCALAR>
dvec<INDEX, SCALAR>::~dvec()
{
  cleanup();
}



// ----------------------------------------------------------------------------
// object management
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
void dvec<INDEX, SCALAR>::resize(INDEX len_)
{
  if (len == len_)
    return;
  
  arraytools::realloc(len_, &X);
  arraytools::check_allocs(X);
  
  if (len_ > len)
    arraytools::zero(len_-len, X+len);
  
  len = len_;
}



template <typename INDEX, typename SCALAR>
void dvec<INDEX, SCALAR>::zero()
{
  if (nnz > 0)
  {
    arraytools::zero(len, X);
    
    nnz = 0;
  }
}



template <typename INDEX, typename SCALAR>
SCALAR dvec<INDEX, SCALAR>::sum() const
{
  SCALAR s = 0;
  #pragma omp simd reduction(+:s)
  for (INDEX i=0; i<len; i++)
    s += X[i];
  
  return s;
}



template <typename INDEX, typename SCALAR>
void dvec<INDEX, SCALAR>::insert(const INDEX i, const SCALAR s)
{
  if (X[i] == 0)
    nnz++;
  else if (s == 0)
    nnz--;
  
  X[i] = s;
}



template <typename INDEX, typename SCALAR>
void dvec<INDEX, SCALAR>::update_nnz()
{
  nnz = 0;
  for (INDEX i=0; i<len; i++)
  {
    if (X[i])
      nnz++;
  }
}



template <typename INDEX, typename SCALAR>
template <typename INDEX_SRC, typename SCALAR_SRC>
void dvec<INDEX, SCALAR>::set(const INDEX_SRC nnz_, const INDEX_SRC *I_, const SCALAR_SRC *X_)
{
  zero();
  
  INDEX_SRC top = I_[nnz_-1];
  if (top > len)
    resize(top + 1);
  
  for (INDEX_SRC i=0; i<nnz_; i++)
    X[I_[i]] = X_[i];
  
  nnz = nnz_;
}

template <typename INDEX, typename SCALAR>
void dvec<INDEX, SCALAR>::set(const spvec<INDEX, SCALAR> &x)
{
  set(x.get_nnz(), x.index_ptr(), x.data_ptr());
}



// ----------------------------------------------------------------------------
// printer
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
void dvec<INDEX, SCALAR>::print() const
{
  printf("## Length %d dense vector with nnz=%d\n", len, nnz);
  for (INDEX i=0; i<len; i++)
    printf("%f\n", (float) X[i]);
}



// ----------------------------------------------------------------------------
// setter/getter
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
const SCALAR dvec<INDEX, SCALAR>::operator[](INDEX i) const
{
  return X[i];
}

template <typename INDEX, typename SCALAR>
SCALAR& dvec<INDEX, SCALAR>::operator[](INDEX i)
{
  return X[i];
}



// ----------------------------------------------------------------------------
// internals
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
void dvec<INDEX, SCALAR>::cleanup()
{
  arraytools::free(X);
  X = NULL;
  
  nnz = 0;
  len = 0;
}


#endif
