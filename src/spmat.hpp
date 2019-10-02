#ifndef SPMAT_CLASS_H
#define SPMAT_CLASS_H


#include <iostream>

#include "arraytools.hpp"
#include "spvec.hpp"


template <typename INDEX, typename SCALAR>
class spmat
{
  public:
    spmat(int nrows_, int ncols_, int len_);
    ~spmat();
    
    void resize(int len_);
    
    void print(bool actual=false) const;
    int insert(const int col, const spvec<INDEX, SCALAR> &x);
    
    int nrows() const {return m;};
    int ncols() const {return n;};
    int get_nnz() const {return nnz;};
    int get_len() const {return len;};
    INDEX* index_ptr() {return I;};
    INDEX* index_ptr() const {return I;};
    INDEX* col_ptr() {return P;};
    INDEX* col_ptr() const {return P;};
    SCALAR* data_ptr() {return X;};
    SCALAR* data_ptr() const {return X;};
  
  protected:
    int m;
    int n;
    int nnz;
    int len;
    int plen;
    INDEX *I;
    INDEX *P;
    SCALAR *X;
  
  private:
    void cleanup();
    void insert_from_ind(const int insertion_ind, const INDEX i, const SCALAR s);
};



template <typename INDEX, typename SCALAR>
spmat<INDEX, SCALAR>::spmat(int nrows_, int ncols_, int len_)
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



template <typename INDEX, typename SCALAR>
void spvec<INDEX, SCALAR>::cleanup()
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


#endif
