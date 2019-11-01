#ifndef DVEC_CLASS_H
#define DVEC_CLASS_H


#include <iostream>

#include "arraytools.hpp"


template <typename INDEX, typename SCALAR>
class dvec
{
  public:
    dvec();
    dvec(INDEX len_);
    ~dvec();
    
    void resize(int len_);
    void zero();
    void insert(const INDEX i, const SCALAR s);
    void update_nnz();
    
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
void dvec<INDEX, SCALAR>::resize(int len_)
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
void dvec<INDEX, SCALAR>::insert(const INDEX i, const SCALAR s)
{
  if (X[i] == 0)
    nnz++;
  
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
