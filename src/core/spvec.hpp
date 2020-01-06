// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_CORE_SPVEC_H
#define SPAR_CORE_SPVEC_H
#pragma once


#include <iostream>

#include "arraytools/src/arraytools.hpp"


template <typename INDEX, typename SCALAR>
class dvec;

template <typename INDEX, typename SCALAR>
class spvec
{
  public:
    spvec();
    spvec(INDEX len_);
    ~spvec();
    
    void resize(INDEX len_);
    void zero();
    INDEX insert(const INDEX i, const SCALAR s);
    template <typename INDEX_SRC, typename SCALAR_SRC>
    void set(const INDEX nnz_, const INDEX_SRC *I_, const SCALAR_SRC *X_);
    void set(const dvec<INDEX, SCALAR> &d);
    
    void print(bool actual=false) const;
    
    INDEX add(const spvec &x);
    INDEX add(const SCALAR *x, const INDEX xlen);
    
    void densify(dvec<INDEX, SCALAR> &d);
    
    INDEX get_nnz() const {return nnz;};
    INDEX get_len() const {return len;};
    INDEX* index_ptr() {return I;};
    INDEX* index_ptr() const {return I;};
    SCALAR* data_ptr() {return X;};
    SCALAR* data_ptr() const {return X;};
  
  protected:
    INDEX nnz;
    INDEX len;
    INDEX *I;
    SCALAR *X;
  
  private:
    void cleanup();
    void insert_from_ind(const INDEX insertion_ind, const INDEX i, const SCALAR s);
};



// ----------------------------------------------------------------------------
// constructor/destructor
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
spvec<INDEX, SCALAR>::spvec()
{
  I = NULL;
  X = NULL;
  
  nnz = 0;
  len = 0;
}



template <typename INDEX, typename SCALAR>
spvec<INDEX, SCALAR>::spvec(INDEX len_)
{
  arraytools::zero_alloc(len_, &I);
  arraytools::zero_alloc(len_, &X);
  
  arraytools::check_allocs(I, X);
  
  nnz = 0;
  len = len_;
}



template <typename INDEX, typename SCALAR>
spvec<INDEX, SCALAR>::~spvec()
{
  cleanup();
}



// ----------------------------------------------------------------------------
// object management
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
void spvec<INDEX, SCALAR>::resize(INDEX len_)
{
  if (len == len_)
    return;
  
  arraytools::realloc(len_, &I);
  arraytools::realloc(len_, &X);
  
  arraytools::check_allocs(I, X);
  
  if (len_ > len)
  {
    arraytools::zero(len_-len, I+len);
    arraytools::zero(len_-len, X+len);
  }
  
  len = len_;
}



template <typename INDEX, typename SCALAR>
void spvec<INDEX, SCALAR>::zero()
{
  if (nnz > 0)
  {
    arraytools::zero(nnz, I);
    arraytools::zero(nnz, X);
    
    nnz = 0;
  }
}



template <typename INDEX, typename SCALAR>
INDEX spvec<INDEX, SCALAR>::insert(const INDEX i, const SCALAR s)
{
  if (nnz == len)
    return 1;
  
  INDEX insertion_ind;
  for (insertion_ind=0; insertion_ind<nnz; insertion_ind++)
  {
    if (i < I[insertion_ind])
      break;
  }
  
  insert_from_ind(insertion_ind, i, s);
  return 0;
}



template <typename INDEX, typename SCALAR>
template <typename INDEX_SRC, typename SCALAR_SRC>
void spvec<INDEX, SCALAR>::set(const INDEX nnz_, const INDEX_SRC *I_, const SCALAR_SRC *X_)
{
  if (len < nnz_)
    resize(nnz_);
  else if (nnz > nnz_)
  {
    arraytools::zero(nnz-nnz_, I+nnz_);
    arraytools::zero(nnz-nnz_, X+nnz_);
  }
  
  arraytools::copy(nnz_, I_, I);
  arraytools::copy(nnz_, X_, X);
  
  nnz = nnz_;
}



template <typename INDEX, typename SCALAR>
void spvec<INDEX, SCALAR>::set(const dvec<INDEX, SCALAR> &d)
{
  INDEX dnnz = d.get_nnz();
  if (dnnz > len)
    resize(dnnz);
  
  INDEX pos = 0;
  SCALAR *d_p = d.data_ptr();
  for (INDEX i=0; i<d.get_len(); i++)
  {
    if (d_p[i])
    {
      I[pos] = i;
      X[pos] = d_p[i];
      
      pos++;
    }
  }
  
  nnz = dnnz;
}



// ----------------------------------------------------------------------------
// printer
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
void spvec<INDEX, SCALAR>::print(bool actual) const
{
  printf("## Length %d sparse vector with nnz=%d\n", len, nnz);
  
  if (actual)
  {
    printf("I: ");
    for (INDEX ind=0; ind<len; ind++)
      std::cout << I[ind] << " ";
    
    printf("\nX: ");
    for (INDEX ind=0; ind<len; ind++)
      std::cout << X[ind] << " ";
    
    putchar('\n');
  }
  else
  {
    const INDEX topval = I[nnz - 1];
    
    INDEX ind = 0;
    for (INDEX ind_g=0; ind_g<=topval; ind_g++)
    {
      if (ind < nnz && ind_g == I[ind])
      {
        std::cout << X[ind] << " ";
        ind++;
      }
      else
        std::cout << (SCALAR) 0 << " ";
    }
    
    std::cout << std::endl;
  }
}



// ----------------------------------------------------------------------------
// adders
// ----------------------------------------------------------------------------

// return needed size of realloc
template <typename INDEX, typename SCALAR>
INDEX spvec<INDEX, SCALAR>::add(const spvec &x)
{
  const INDEX *xI = x.index_ptr();
  const SCALAR *xX = x.data_ptr();
  INDEX ind = 0;
  
  // pre-scan to see if a re-alloc is necessary
  INDEX num_inserted = 0;
  for (INDEX xind=0; xind<x.get_nnz(); xind++)
  {
    const INDEX xi = xI[xind];
    while (ind < nnz && I[ind] < xi)
      ind++;
    
    if (ind >= nnz || I[ind] > xi)
      num_inserted++;
  }
  
  if (num_inserted > (len - nnz))
    return num_inserted;
  
  // add the vectors
  ind = 0;
  for (INDEX xind=0; xind<x.get_nnz(); xind++)
  {
    const INDEX xi = xI[xind];
    while (ind < nnz && I[ind] < xi)
      ind++;
    
    if (I[ind] == xi)
      X[ind++] += xX[xind];
    else if (ind == nnz || I[ind] > xi)
      insert_from_ind(ind, xi, xX[xind]);
  }
  
  return 0;
}



template <typename INDEX, typename SCALAR>
INDEX spvec<INDEX, SCALAR>::add(const SCALAR *x, const INDEX xlen)
{
  INDEX ind = 0;
  
  // pre-scan to see if a re-alloc is necessary
  INDEX num_inserted = 0;
  for (INDEX xi=0; xi<xlen; xi++)
  {
    if (x[xi] == (SCALAR)0)
      continue;
    
    while (ind < nnz && I[ind] < xi)
      ind++;
    
    if (ind >= nnz || I[ind] > xi)
      num_inserted++;
  }
  
  if (num_inserted > (xlen - nnz))
    return num_inserted;
  
  // add the vectors
  ind = 0;
  for (INDEX xi=0; xi<xlen; xi++)
  {
    if (x[xi] == (SCALAR)0)
      continue;
    
    while (ind < nnz && I[ind] < xi)
      ind++;
    
    if (I[ind] == xi)
      X[ind++] += x[xi];
    else if (ind == nnz || I[ind] > xi)
      insert_from_ind(ind, xi, x[xi]);
  }
  
  return 0;
}



// ----------------------------------------------------------------------------
// converters
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
void spvec<INDEX, SCALAR>::densify(dvec<INDEX, SCALAR> &d)
{
  if (I[nnz-1] > d.get_len())
    throw std::logic_error("dense array not large enough to store sparse vector");
  
  d.zero();
  for (INDEX pos=0; pos<nnz; pos++)
    d.insert(I[pos], X[pos]);
}



// ----------------------------------------------------------------------------
// internals
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
void spvec<INDEX, SCALAR>::cleanup()
{
  arraytools::free(I);
  I = NULL;
  
  arraytools::free(X);
  X = NULL;
  
  nnz = 0;
  len = 0;
}



template <typename INDEX, typename SCALAR>
void spvec<INDEX, SCALAR>::insert_from_ind(const INDEX insertion_ind,
  const INDEX i, const SCALAR s)
{
  for (INDEX ind=nnz; ind>insertion_ind; ind--)
  {
    I[ind] = I[ind-1];
    X[ind] = X[ind-1];
  }
  
  I[insertion_ind] = i;
  X[insertion_ind] = s;
  
  nnz++;
}


#endif
