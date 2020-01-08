#ifndef SPAR_TESTS_MPI_GEN_H
#define SPAR_TESTS_MPI_GEN_H
#pragma once


template <typename INDEX, typename SCALAR>
static inline void fill_sparse_mat(spmat<INDEX, SCALAR> &x)
{
  spvec<INDEX, SCALAR> s(3);
  
  s.zero();
  s.insert(0, 1);
  s.insert(5, 1);
  s.insert(9, 1);
  x.insert(0, s);
  
  s.zero();
  s.insert(1, 2);
  s.insert(3, 1);
  x.insert(2, s);
  
  s.zero();
  s.insert(2, 2);
  s.insert(4, 1);
  x.insert(6, s);
  
  if (rank != 0)
  {
    s.zero();
    s.insert(5, 1);
    x.insert(5, s);
  }
}


#endif
