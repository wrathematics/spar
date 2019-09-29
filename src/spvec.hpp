#ifndef SPVEC_CLASS_H
#define SPVEC_CLASS_H


#include <iostream>

#include "arraytools.hpp"


template <typename INDEX, typename SCALAR>
class spvec
{
  public:
    spvec(int len_);
    ~spvec();
    
    void resize(int len_);
    void set(int nnz_, INDEX *I_, SCALAR *X_);
    void zero();
    
    void print(bool actual=false) const;
    int insert(const INDEX i, const SCALAR s);
    int add(const spvec &x);
    int add(const SCALAR *x, const int xlen);
    
    int get_nnz() const {return nnz;};
    int get_len() const {return len;};
    INDEX* index_ptr() {return I;};
    INDEX* index_ptr() const {return I;};
    SCALAR* data_ptr() {return X;};
    SCALAR* data_ptr() const {return X;};
  
  protected:
    int nnz;
    int len;
    INDEX *I;
    SCALAR *X;
  
  private:
    void cleanup();
    void insert_from_ind(const int insertion_ind, const INDEX i, const SCALAR s);
};



template <typename INDEX, typename SCALAR>
spvec<INDEX, SCALAR>::spvec(int len_)
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



template <typename INDEX, typename SCALAR>
void spvec<INDEX, SCALAR>::resize(int len_)
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
void spvec<INDEX, SCALAR>::set(int nnz_, INDEX *I_, SCALAR *X_)
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
void spvec<INDEX, SCALAR>::print(bool actual) const
{
  printf("## Length %d sparse vector with nnz=%d\n", len, nnz);
  
  if (actual)
  {
    printf("I: ");
    for (int ind=0; ind<len; ind++)
      std::cout << I[ind] << " ";
    
    printf("\nX: ");
    for (int ind=0; ind<len; ind++)
      std::cout << X[ind] << " ";
    
    putchar('\n');
  }
  else
  {
    int ind_s = 0;
    for (int ind=0; ind<len; ind++)
    {
      if (ind_s < nnz && ind == I[ind_s])
      {
        std::cout << X[ind_s];
        ind_s++;
      }
      else
        std::cout << (SCALAR) 0;
      
      std::cout << " ";
    }
    
    std::cout << std::endl;
  }
}



template <typename INDEX, typename SCALAR>
int spvec<INDEX, SCALAR>::insert(const INDEX i, const SCALAR s)
{
  if (nnz == len)
    return 1;
  
  int insertion_ind;
  for (insertion_ind=0; insertion_ind<nnz; insertion_ind++)
  {
    if (i < I[insertion_ind])
      break;
  }
  
  insert_from_ind(insertion_ind, i, s);
  return 0;
}



// return needed size of realloc
template <typename INDEX, typename SCALAR>
int spvec<INDEX, SCALAR>::add(const spvec &x)
{
  const INDEX *xI = x.index_ptr();
  const SCALAR *xX = x.data_ptr();
  int ind = 0;
  
  // pre-scan to see if a re-alloc is necessary
  int num_inserted = 0;
  for (int xind=0; xind<x.get_nnz(); xind++)
  {
    const int xi = xI[xind];
    while (I[ind] < xi)
      ind++;
    
    if (I[ind] > xi)
      num_inserted++;
  }
  
  if (num_inserted > (len - nnz))
    return num_inserted;
  
  // add the vectors
  ind = 0;
  for (int xind=0; xind<x.get_nnz(); xind++)
  {
    const int xi = xI[xind];
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
int spvec<INDEX, SCALAR>::add(const SCALAR *x, const int xlen)
{
  int ind = 0;
  
  // pre-scan to see if a re-alloc is necessary
  int num_inserted = 0;
  for (int xi=0; xi<xlen; xi++)
  {
    if (x[xi] == (SCALAR)0)
      continue;
    
    while (I[ind] < xi)
      ind++;
    
    if (I[ind] > xi)
      num_inserted++;
  }
  
  if (num_inserted > (xlen - nnz))
    return num_inserted;
  
  // add the vectors
  ind = 0;
  for (int xi=0; xi<xlen; xi++)
  {
    while (ind < nnz && I[ind] < xi)
      ind++;
    
    if (I[ind] == xi)
      X[ind++] += x[xi];
    else if (ind == nnz || I[ind] > xi)
      insert_from_ind(ind, xi, x[xi]);
  }
  
  return 0;
}



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
void spvec<INDEX, SCALAR>::insert_from_ind(const int insertion_ind,
  const INDEX i, const SCALAR s)
{
  for (int ind=nnz; ind>insertion_ind; ind--)
  {
    I[ind] = I[ind-1];
    X[ind] = X[ind-1];
  }
  
  I[insertion_ind] = i;
  X[insertion_ind] = s;
  
  nnz++;
}


#endif
