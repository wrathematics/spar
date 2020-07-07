// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_CORE_DVEC_H
#define SPAR_CORE_DVEC_H
#pragma once


#include <iostream>

#include "../arraytools/src/arraytools.hpp"


namespace spar
{
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
      void insert(const INDEX i, const SCALAR s);
      void update_nnz();
      
      void print() const;
      void info() const;
      
      SCALAR sum() const;
      const SCALAR operator[](INDEX i) const; // getter
      SCALAR& operator[](INDEX i); // setter
      
      template <typename INDEX_SRC, typename SCALAR_SRC>
      void set(const INDEX_SRC nnz_, const INDEX_SRC *I_, const SCALAR_SRC *X_);
      void set(const spvec<INDEX, SCALAR> &x);
      
      /// Number of non-zero elements.
      INDEX get_nnz() const {return nnz;};
      /// Length of the index and data arrays.
      INDEX get_len() const {return len;};
      /// Return a pointer to the data array `X`.
      SCALAR* data_ptr() {return X;};
      /// \overload
      SCALAR* data_ptr() const {return X;};
    
    protected:
      INDEX nnz;
      INDEX len;
      SCALAR *X;
    
    private:
      void cleanup();
  };
}



// ----------------------------------------------------------------------------
// constructor/destructor
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
spar::dvec<INDEX, SCALAR>::dvec()
{
  X = NULL;
  
  nnz = 0;
  len = 0;
}



/**
  @brief Constructor.
  
  @param[in] len_ Length of the vector (elements, not bytes).
  
  @allocs One internal array is allocated.
  
  @except If a memory allocation fails, a `bad_alloc` exception will be thrown.
 */
template <typename INDEX, typename SCALAR>
spar::dvec<INDEX, SCALAR>::dvec(INDEX len_)
{
  arraytools::zero_alloc(len_, &X);
  arraytools::check_allocs(X);
  
  nnz = 0;
  len = len_;
}



template <typename INDEX, typename SCALAR>
spar::dvec<INDEX, SCALAR>::~dvec()
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
void spar::dvec<INDEX, SCALAR>::resize(INDEX len_)
{
  if (len == len_)
    return;
  
  arraytools::realloc(len_, &X);
  arraytools::check_allocs(X);
  
  if (len_ > len)
    arraytools::zero(len_-len, X+len);
  
  len = len_;
}



/// Zero all data in the dense vector. Performs no allocations or resizing.
template <typename INDEX, typename SCALAR>
void spar::dvec<INDEX, SCALAR>::zero()
{
  if (nnz > 0)
  {
    arraytools::zero(len, X);
    
    nnz = 0;
  }
}



/**
  @brief Insert a value at the specified index.
  
  @param[in] i The index.
  @param[in] s The input value.
 */
template <typename INDEX, typename SCALAR>
void spar::dvec<INDEX, SCALAR>::insert(const INDEX i, const SCALAR s)
{
  if (X[i] == 0)
    nnz++;
  else if (s == 0)
    nnz--;
  
  X[i] = s;
}



/**
  @brief Updates the internal "number non-zero" count. Useful if operating
  directly on the array.
 */
template <typename INDEX, typename SCALAR>
void spar::dvec<INDEX, SCALAR>::update_nnz()
{
  nnz = 0;
  for (INDEX i=0; i<len; i++)
  {
    if (X[i])
      nnz++;
  }
}



// ----------------------------------------------------------------------------
// printer
// ----------------------------------------------------------------------------

/// Print the vector.
template <typename INDEX, typename SCALAR>
void spar::dvec<INDEX, SCALAR>::print() const
{
  for (INDEX i=0; i<len; i++)
    printf("%f\n", (float) X[i]);
}



/// Print some quick info about the sparse vector.
template <typename INDEX, typename SCALAR>
void spar::dvec<INDEX, SCALAR>::info() const
{
  printf("# dvec");
  printf(" %d", len);
  printf(" with nnz=%d", nnz);
  printf(" (index=%s scalar=%s)", typeid(INDEX).name(), typeid(SCALAR).name());
  printf("\n");
}



// ----------------------------------------------------------------------------
// setter/getter
// ----------------------------------------------------------------------------

/**
  @brief Sum the vector elements.
  
  @return The sum of the vector elements.
 */
template <typename INDEX, typename SCALAR>
SCALAR spar::dvec<INDEX, SCALAR>::sum() const
{
  SCALAR s = 0;
  #pragma omp simd reduction(+:s)
  for (INDEX i=0; i<len; i++)
    s += X[i];
  
  return s;
}



/**
  @brief Getter.
  
  @param[in] i Desired index.
  
  @return Returns the scalar at index `i`.
 */
template <typename INDEX, typename SCALAR>
const SCALAR spar::dvec<INDEX, SCALAR>::operator[](INDEX i) const
{
  return X[i];
}



/**
  @brief Setter. The RHS is the scalar value to set the vector at index `i`.
  
  @param[in] i Desired index.
 */
template <typename INDEX, typename SCALAR>
SCALAR& spar::dvec<INDEX, SCALAR>::operator[](INDEX i)
{
  return X[i];
}



// ----------------------------------------------------------------------------
// converters
// ----------------------------------------------------------------------------

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
void spar::dvec<INDEX, SCALAR>::set(const INDEX_SRC nnz_, const INDEX_SRC *I_, const SCALAR_SRC *X_)
{
  zero();
  
  INDEX_SRC top = I_[nnz_-1];
  if (top > len)
    resize(top + 1);
  
  for (INDEX_SRC i=0; i<nnz_; i++)
    X[I_[i]] = X_[i];
  
  nnz = nnz_;
}



/// \overload
template <typename INDEX, typename SCALAR>
void spar::dvec<INDEX, SCALAR>::set(const spvec<INDEX, SCALAR> &x)
{
  set(x.get_nnz(), x.index_ptr(), x.data_ptr());
}



// ----------------------------------------------------------------------------
// internals
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
void spar::dvec<INDEX, SCALAR>::cleanup()
{
  arraytools::free(X);
  X = NULL;
  
  nnz = 0;
  len = 0;
}


#endif
