// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_CORE_SPMAT_H
#define SPAR_CORE_SPMAT_H
#pragma once


#include <iostream>

#include "../arraytools/src/arraytools.hpp"
#include "defs.hpp"


template <typename INDEX, typename SCALAR>
class spvec;

/**
  @brief Basic sparse matrix class in CSC format.
  
  @tparam INDEX should be some kind of fundamental indexing type, like `int`
  or `uint16_t`.
  @tparam SCALAR should be a fundamental numeric type like `int` or `float`.
 */
template <typename INDEX, typename SCALAR>
class spmat
{
  public:
    spmat(INDEX nrows_, INDEX ncols_, INDEX len_);
    ~spmat();
    
    void resize(INDEX len_);
    void zero();
    INDEX insertable(const spvec<INDEX, SCALAR> &x);
    void insert(const INDEX col, const spvec<INDEX, SCALAR> &x);
    void update_nnz();
    void get_col(const INDEX col, spvec<INDEX, SCALAR> &x) const;
    
    void print(bool actual=false);
    void info() const;
    
    float sparsity() const;
    float density() const;
    
    /// Number of rows.
    INDEX nrows() const {return m;};
    /// Number of columns.
    INDEX ncols() const {return n;};
    /// Number of non-zero elements.
    INDEX get_nnz() const {return nnz;};
    /// Length of the index and data arrays.
    INDEX get_len() const {return len;};
    /// Return a pointer to the index array `I`.
    INDEX* index_ptr() {return I;};
    /// \overload
    INDEX* index_ptr() const {return I;};
    /// Return a pointer to the column array `P`.
    INDEX* col_ptr() {return P;};
    /// \overload
    INDEX* col_ptr() const {return P;};
    /// Return a pointer to the data array `X`.
    SCALAR* data_ptr() {return X;};
    /// \overload
    SCALAR* data_ptr() const {return X;};
  
  protected:
    INDEX m;
    INDEX n;
    INDEX nnz;
    INDEX len;
    INDEX plen;
    INDEX *I;
    INDEX *P;
    SCALAR *X;
  
  private:
    void cleanup();
    INDEX* csc2coo();
    void insert_spvec(const INDEX col, const spvec<INDEX, SCALAR> &x);
};



// ----------------------------------------------------------------------------
// constructor/destructor
// ----------------------------------------------------------------------------

/**
  @brief Constructor.
  
  @param[in] nrows_,ncols_ The dimension of the matrix.
  @param[in] len_ The amount of storage to initially allocate (elements, not
  bytes).
  
  @allocs Three internal arrays are allocated.
  
  @except If a memory allocation fails, a `bad_alloc` exception will be thrown.
 */
template <typename INDEX, typename SCALAR>
spmat<INDEX, SCALAR>::spmat(INDEX nrows_, INDEX ncols_, INDEX len_)
{
  arraytools::zero_alloc(len_, &I);
  arraytools::zero_alloc(ncols_+1, &P);
  arraytools::zero_alloc(len_, &X);
  
  arraytools::check_allocs(I, P, X);
  
  m = nrows_;
  n = ncols_;
  
  nnz = 0;
  len = len_;
  plen = ncols_ + 1;
}



template <typename INDEX, typename SCALAR>
spmat<INDEX, SCALAR>::~spmat()
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
void spmat<INDEX, SCALAR>::resize(INDEX len_)
{
  if (len == len_)
    return;
  
  arraytools::realloc(len_, &I);
  arraytools::realloc(len_, &X);
  
  arraytools::check_allocs(I, P, X);
  
  len = len_;
}



/// Zero all data in the sparse matrix. Performs no allocations or resizing.
template <typename INDEX, typename SCALAR>
void spmat<INDEX, SCALAR>::zero()
{
  if (nnz > 0)
  {
    arraytools::zero(len, I);
    arraytools::zero(n+1, P);
    arraytools::zero(len, X);
    
    nnz = 0;
  }
}



/**
  @brief Check if the matrix has enough storage capacity to hold the given
  sparse vector.
  
  @param[in] x The input column.
  
  @return The number of elements needed that exceed the current matrix storage
  capacity (0 if the vector fits).
 */
template <typename INDEX, typename SCALAR>
INDEX spmat<INDEX, SCALAR>::insertable(const spvec<INDEX, SCALAR> &x)
{
  if (x.get_nnz() > len - nnz)
    return x.get_nnz() - (len - nnz);
  else
    return (INDEX) 0;
}



/**
  @brief Insert a sparse vector into the specified column.
  
  @param[in] col The column index.
  @param[in] x The input column.
  
  @allocs The internal arrays will resize themselves as needed.
  
  @except If a memory allocation fails, a `bad_alloc` exception will be thrown.
 */
template <typename INDEX, typename SCALAR>
void spmat<INDEX, SCALAR>::insert(const INDEX col, const spvec<INDEX, SCALAR> &x)
{
  INDEX needed_space = len - nnz;
  if (x.get_nnz() > needed_space)
    resize((len + needed_space) * spar::internal::defs::MEM_FUDGE_ELT_FAC);
  
  insert_spvec(col, x);
}



/**
  @brief Updates the internal "number non-zero" count. Useful if operating
  directly on the internal arrays.
 */
template <typename INDEX, typename SCALAR>
void spmat<INDEX, SCALAR>::update_nnz()
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
  
  @param[in] col The column index.
  @param[out] x The input column.
  
  @allocs The return vector `x` will be resized as needed.
  
  @except If a memory allocation fails, a `bad_alloc` exception will be thrown.
 */
template <typename INDEX, typename SCALAR>
void spmat<INDEX, SCALAR>::get_col(const INDEX col, spvec<INDEX, SCALAR> &x) const
{
  const INDEX ind = P[col];
  if (P[col + 1] == ind)
  {
    x.zero();
    return;
  }
  
  const INDEX col_nnz = P[col + 1] - ind;
  x.set(col_nnz, I + ind, X + ind);
}



// ----------------------------------------------------------------------------
// utils
// ----------------------------------------------------------------------------

/**
  @brief Return the proportion of sparsity.
  
  @return The non-zero elements divided by the matrix dimensions.
 */
template <typename INDEX, typename SCALAR>
float spmat<INDEX, SCALAR>::sparsity() const
{
  return (float)nnz/m/n;
}



/**
  @brief Return the proportion of density.
  
  @return The complement of the number of non-zero elements divided by the
  matrix dimensions from one (i.e., 1 minus the sparsity).
 */
template <typename INDEX, typename SCALAR>
float spmat<INDEX, SCALAR>::density() const
{
  return 1.f - sparsity();
}



// ----------------------------------------------------------------------------
// printer
// ----------------------------------------------------------------------------

/**
  @brief Print the matrix.
  
  @param[in] actual Should the actual/literal storage (some arrays) be printed?
  Otherwise, the conceptual dense matrix will be printed.
  
  @allocs If not printing the actual values, a temporary storage array is needed
  so as to convert the matrix from CSC to COO. It is of length
  `len*sizeof(INDEX)`.
  
  @except If a memory allocation fails, a `bad_alloc` exception will be thrown.
 */
template <typename INDEX, typename SCALAR> // NOTE to self: don't set this method const
void spmat<INDEX, SCALAR>::print(bool actual)
{
  if (actual)
  {
    printf("I: ");
    for (INDEX ind=0; ind<len; ind++)
      std::cout << I[ind] << " ";
    
    printf("\nP: ");
    for (INDEX ind=0; ind<plen; ind++)
      std::cout << P[ind] << " ";
    
    printf("\nX: ");
    for (INDEX ind=0; ind<len; ind++)
      std::cout << X[ind] << " ";
    
    putchar('\n');
  }
  else
  {
    if (nnz == 0)
    {
      for (INDEX i=0; i<m; i++)
      {
        for (INDEX j=0; j<n; j++)
          std::cout << (SCALAR) 0 << " ";
        
        putchar('\n');
      }
    }
    else
    {
      INDEX* CI = csc2coo();
      for (INDEX i=0; i<m; i++)
      {
        INDEX ind = 0;
        for (INDEX j=0; j<n; j++)
        {
          while (I[ind] < i || CI[ind] < j)
            ind++;
          
          if (I[ind] == i && CI[ind] == j)
            std::cout << X[ind++];
          else
            std::cout << (SCALAR) 0;
          
          putchar(' ');
        }
        
        putchar('\n');
      }
      
      arraytools::free(CI);
    }
  }
  
  putchar('\n');
}



/// Print some quick info about the sparse matrix.
template <typename INDEX, typename SCALAR>
void spmat<INDEX, SCALAR>::info() const
{
  printf("# spmat");
  printf(" %dx%d", m, n);
  printf(" with nnz=%d (%.2f%% sparse)", nnz, sparsity()*100.f);
  printf(" and len=%d", len);
  printf(" (index=%s scalar=%s)", typeid(INDEX).name(), typeid(SCALAR).name());
  printf("\n");
}



// ----------------------------------------------------------------------------
// internals
// ----------------------------------------------------------------------------

template <typename INDEX, typename SCALAR>
void spmat<INDEX, SCALAR>::cleanup()
{
  arraytools::free(I);
  I = NULL;
  
  arraytools::free(P);
  P = NULL;
  
  arraytools::free(X);
  X = NULL;
  
  m = 0;
  n = 0;
  
  nnz = 0;
  len = 0;
}



// convert "column pointer" P into column index
template <typename INDEX, typename SCALAR>
INDEX* spmat<INDEX, SCALAR>::csc2coo()
{
  INDEX j = 0;
  INDEX ind = 0;
  
  INDEX* CI;
  arraytools::alloc(len, &CI);
  if (CI == NULL)
  {
    cleanup();
    throw std::bad_alloc();
  }
  
  for (INDEX c=0; c<plen; c++)
  {
    INDEX diff = P[c+1] - P[c];
    
    while (diff > 0)
    {
      CI[ind] = j;
      
      ind++;
      diff--;
    }
    
    j++;
  }
  
  return CI;
}



template <typename INDEX, typename SCALAR>
void spmat<INDEX, SCALAR>::insert_spvec(const INDEX col, const spvec<INDEX, SCALAR> &x)
{
  const INDEX *xI = x.index_ptr();
  const SCALAR *xX = x.data_ptr();
  
  INDEX ind = P[col];
  const INDEX xnnz = x.get_nnz();
  for (INDEX xind=0; xind<xnnz; xind++)
  {
    I[ind] = xI[xind];
    X[ind] = xX[xind];
    
    ind++;
  }
  
  for (INDEX ind=col+1; ind<plen; ind++)
    P[ind] += xnnz;
  
  nnz += xnnz;
}


#endif
