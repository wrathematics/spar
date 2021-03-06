// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_CORE_SPVEC_H
#define SPAR_CORE_SPVEC_H
#pragma once


#include <iostream>

#include "../arraytools/src/arraytools.hpp"


namespace spar
{
  template <typename INDEX, typename SCALAR>
  class dvec;
  
  /**
    @brief Basic sparse vector class.
    
    @tparam INDEX should be some kind of fundamental indexing type, like `int`
    or `uint16_t`.
    @tparam SCALAR should be a fundamental numeric type like `int` or `float`.
   */
  template <typename INDEX, typename SCALAR>
  class spvec
  {
    public:
      spvec();
      spvec(INDEX len_);
      ~spvec();
      
      void resize(INDEX len_);
      void zero();
      INDEX insertable() const;
      void insert(const INDEX i, const SCALAR s);
      void update_nnz();
      SCALAR get(const INDEX ind) const;
      
      void print(bool actual=false) const;
      void info() const;
      
      INDEX add(const spvec &x);
      INDEX add(const SCALAR *x, const INDEX xlen);
      
      void densify(dvec<INDEX, SCALAR> &d) const;
      template <typename INDEX_SRC, typename SCALAR_SRC>
      void set(const INDEX nnz_, const INDEX_SRC *I_, const SCALAR_SRC *X_);
      void set(const dvec<INDEX, SCALAR> &d);
      
      /// Number of non-zero elements.
      INDEX get_nnz() const {return nnz;};
      /// Length of the index and data arrays.
      INDEX get_len() const {return len;};
      /// Return a pointer to the index array `I`.
      INDEX* index_ptr() {return I;};
      /// \overload
      INDEX* index_ptr() const {return I;};
      /// Return a pointer to the data array `X`.
      SCALAR* data_ptr() {return X;};
      /// \overload
      SCALAR* data_ptr() const {return X;};
    
    protected:
      /// Number non-zero.
      INDEX nnz;
      /// Length of internal row/data arrays.
      INDEX len;
      /// Index array.
      INDEX *I;
      /// Data array.
      SCALAR *X;
    
    private:
      void cleanup();
      void insert_from_ind(const INDEX insertion_ind, const INDEX i, const SCALAR s);
  };
}



// ----------------------------------------------------------------------------
// constructor/destructor
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
spar::spvec<INDEX, SCALAR>::spvec()
{
  I = NULL;
  X = NULL;
  
  nnz = 0;
  len = 0;
}



/**
  @brief Constructor.
  
  @param[in] len_ The amount of storage to initially allocate (elements, not
  bytes).
  
  @allocs Two internal arrays are allocated.
  
  @except If a memory allocation fails, a `bad_alloc` exception will be thrown.
 */
template <typename INDEX, typename SCALAR>
spar::spvec<INDEX, SCALAR>::spvec(INDEX len_)
{
  arraytools::zero_alloc(len_, &I);
  arraytools::zero_alloc(len_, &X);
  
  arraytools::check_allocs(I, X);
  
  nnz = 0;
  len = len_;
}



template <typename INDEX, typename SCALAR>
spar::spvec<INDEX, SCALAR>::~spvec()
{
  cleanup();
}



// ----------------------------------------------------------------------------
// object management
// ----------------------------------------------------------------------------

/**
  @brief Resize the internal storage.
  
  @param[in] len_ The new amount of internal storage to use (elements, not
  bytes).
  
  @allocs The internal arrays will resize themselves as needed.
  
  @except If a memory allocation fails, a `bad_alloc` exception will be thrown.
 */
template <typename INDEX, typename SCALAR>
void spar::spvec<INDEX, SCALAR>::resize(INDEX len_)
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



/// Zero all data in the sparse vector. Performs no allocations or resizing.
template <typename INDEX, typename SCALAR>
void spar::spvec<INDEX, SCALAR>::zero()
{
  if (nnz > 0)
  {
    arraytools::zero(nnz, I);
    arraytools::zero(nnz, X);
    
    nnz = 0;
  }
}



/**
  @brief Check if the vector has enough storage capacity to hold an additional
  scalar.
  
  @return The number of additional elements the internal vector storage needs
  to store an additional scalar: 0 or 1.
 */
template <typename INDEX, typename SCALAR>
INDEX spar::spvec<INDEX, SCALAR>::insertable() const
{
  if (nnz == len)
    return (INDEX) 1;
  else
    return (INDEX) 0;
}



/**
  @brief Insert a value at the specified index.
  
  @param[in] i The index.
  @param[in] s The input value.
  
  @allocs The internal arrays will resize themselves as needed.
  
  @except If a memory allocation fails, a `bad_alloc` exception will be thrown.
 */
template <typename INDEX, typename SCALAR>
void spar::spvec<INDEX, SCALAR>::insert(const INDEX i, const SCALAR s)
{
  if (nnz == len)
    resize(len + 1);
  
  INDEX insertion_ind;
  for (insertion_ind=0; insertion_ind<nnz; insertion_ind++)
  {
    if (i < I[insertion_ind])
      break;
  }
  
  insert_from_ind(insertion_ind, i, s);
}



/**
  @brief Updates the internal "number non-zero" count. Useful if operating
  directly on the internal arrays.
 */
template <typename INDEX, typename SCALAR>
void spar::spvec<INDEX, SCALAR>::update_nnz()
{
  nnz = 0;
  for (INDEX i=0; i<len; i++)
  {
    if (X[i])
      nnz++;
    else
      break;
  }
}



/**
  @brief Retrieve the specified column as a sparse vector.
  
  @param[in] ind The index.
  @return The scalar value at position `ind`. Will return 0 if no non-zero value
  is stored at that index.
 */
template <typename INDEX, typename SCALAR>
SCALAR spar::spvec<INDEX, SCALAR>::get(const INDEX ind) const
{
  SCALAR val = 0;
  
  for (INDEX pos=0; pos<nnz; pos++)
  {
    if (I[pos] == ind)
    {
      val = X[pos];
      break;
    }
  }
  
  return val;
}



// ----------------------------------------------------------------------------
// printer
// ----------------------------------------------------------------------------

/**
  @brief Print the vector.
  
  @param[in] actual Should the actual/literal storage (some arrays) be printed?
  Otherwise, the conceptual dense vector will be printed.
 */
template <typename INDEX, typename SCALAR>
void spar::spvec<INDEX, SCALAR>::print(bool actual) const
{
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



/// Print some quick info about the sparse vector.
template <typename INDEX, typename SCALAR>
void spar::spvec<INDEX, SCALAR>::info() const
{
  printf("# spvec");
  printf(" %d", len);
  printf(" with nnz=%d", nnz);
  printf(" (index=%s scalar=%s)", typeid(INDEX).name(), typeid(SCALAR).name());
  printf("\n");
}



// ----------------------------------------------------------------------------
// adders
// ----------------------------------------------------------------------------

/**
  @brief Add the input to the sparse vector.
  
  @param[in] x Input to adder.
  
  @return Returns the needed size of realloc.
 */
template <typename INDEX, typename SCALAR>
INDEX spar::spvec<INDEX, SCALAR>::add(const spvec &x)
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



/// \overload
template <typename INDEX, typename SCALAR>
INDEX spar::spvec<INDEX, SCALAR>::add(const SCALAR *x, const INDEX xlen)
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

/**
  @brief Convert the sparse vector into a dense vector.
  
  @param[in] d The output dense array.
  
  @allocs The output (dense vector) will be resized as needed.
  
  @except If a memory allocation fails, a `bad_alloc` exception will be thrown.
 */
template <typename INDEX, typename SCALAR>
void spar::spvec<INDEX, SCALAR>::densify(dvec<INDEX, SCALAR> &d) const
{
  if (nnz && I[nnz-1] > d.get_len())
    throw std::logic_error("dense array not large enough to store sparse vector");
  
  d.zero();
  for (INDEX pos=0; pos<nnz; pos++)
    d.insert(I[pos], X[pos]);
}



/**
  @brief Set the vector to the values in the input.
  
  @param[in] nnz_ Length of the input index/scalar arrays.
  @param[in] I_ Index array.
  @param[in] X_ Scalar array.
  
  @allocs The internal arrays will resize themselves as needed.
  
  @except If a memory allocation fails, a `bad_alloc` exception will be thrown.
 */
template <typename INDEX, typename SCALAR>
template <typename INDEX_SRC, typename SCALAR_SRC>
void spar::spvec<INDEX, SCALAR>::set(const INDEX nnz_, const INDEX_SRC *I_, const SCALAR_SRC *X_)
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



/// \overload
template <typename INDEX, typename SCALAR>
void spar::spvec<INDEX, SCALAR>::set(const dvec<INDEX, SCALAR> &d)
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
// internals
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
void spar::spvec<INDEX, SCALAR>::cleanup()
{
  arraytools::free(I);
  I = NULL;
  
  arraytools::free(X);
  X = NULL;
  
  nnz = 0;
  len = 0;
}



template <typename INDEX, typename SCALAR>
void spar::spvec<INDEX, SCALAR>::insert_from_ind(const INDEX insertion_ind,
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
