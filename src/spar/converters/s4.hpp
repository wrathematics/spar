// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_CONVERTERS_S4_H
#define SPAR_CONVERTERS_S4_H
#pragma once


#include <Rdefines.h>
#include <Rinternals.h>

#undef nrows
#undef ncols


namespace spar
{
  template <typename INDEX, typename SCALAR>
  class spvec;

  template <typename INDEX, typename SCALAR>
  class spmat;
  
  namespace internal
  {
    namespace sexp
    {
      static inline SEXP get_obj_from_s4(SEXP s4, const char *obj)
      {
        return GET_SLOT(s4, install(obj));
      }
      
      static inline SEXP get_x_from_s4(SEXP s4)
      {
        return get_obj_from_s4(s4, "x");
      }
      
      static inline SEXP get_i_from_s4(SEXP s4)
      {
        return get_obj_from_s4(s4, "i");
      }
      
      static inline SEXP get_p_from_s4(SEXP s4)
      {
        return get_obj_from_s4(s4, "p");
      }
      
      static inline void get_dim_from_s4(SEXP s4, int *m, int *n)
      {
        SEXP dim = get_obj_from_s4(s4, "Dim");
        *m = INTEGER(dim)[0];
        *n = INTEGER(dim)[1];
      }
      
      static inline int get_nnz_from_s4(SEXP s4_I)
      {
        return LENGTH(s4_I);
      }
      
      static inline int get_plen_from_s4(SEXP s4_P)
      {
        return LENGTH(s4_P);
      }
      
      static inline int get_col_len_from_s4(int col_ind, SEXP s4_P)
      {
        return INTEGER(s4_P)[col_ind+1] - INTEGER(s4_P)[col_ind];
      }
    }
  }
  
  
  
  /// Converters.
  namespace conv
  {
    /**
      @brief Convert a column of a `dgCMatrix` object into an `spvec` sparse
      vector.
      
      @param[in] col_ind Zero-based index of the desired column.
      @param[in] s4 The input `dgCMatrix` object.
      @param[out] s The return sparse matrix.
      
      @allocs The passed sparse vector `s` will resize itself as needed during
      function execution.
      
      @except If a memory allocation fails, a `bad_alloc` exception will be
      thrown.
      
      @tparam INDEX should be some kind of fundamental indexing type, like `int`
      or `uint16_t`.
      @tparam SCALAR should be a fundamental numeric type like `int` or `float`.
     */
    template <typename INDEX, typename SCALAR>
    static inline void s4col_to_spvec(const int col_ind, SEXP s4, spvec<INDEX, SCALAR> &s)
    {
      SEXP s4_X = internal::sexp::get_x_from_s4(s4);
      SEXP s4_I = internal::sexp::get_i_from_s4(s4); // len == nnz
      SEXP s4_P = internal::sexp::get_p_from_s4(s4);
      
      const int start_ind = INTEGER(s4_P)[col_ind];
      const int col_len = INTEGER(s4_P)[col_ind+1] - start_ind;
      
      s.set(col_len, INTEGER(s4_I) + start_ind, REAL(s4_X) + start_ind);
    }
    
    /// \overload
    template <typename INDEX, typename SCALAR>
    static inline spvec<INDEX, SCALAR> s4col_to_spvec(const int col_ind, SEXP s4)
    {
      SEXP s4_P = internal::sexp::get_p_from_s4(s4);
      const int col_len = internal::sexp::get_col_len_from_s4(col_ind, s4_P);
      
      spvec<INDEX, SCALAR> s(col_len * spar::internal::defs::MEM_FUDGE_ELT_FAC);
      s4col_to_spvec(col_ind, s4, s);
      return s;
    }
    
    
    
    /**
      @brief Convert an `spmat` object into a `dgCMatrix` sparse matrix.
      
      @param[in] s The input `spmat` object.
      
      @return The return sparse matrix.
      
      @allocs The return object is roughly of size:
      `sizeof(int)*(2 + nnz + (n+1)) + sizeof(double)*nnz`.
      
      @except If a memory allocation fails, a `bad_alloc` exception will be
      thrown.
      
      @tparam INDEX should be some kind of fundamental indexing type, like `int`
      or `uint16_t`.
      @tparam SCALAR should be a fundamental numeric type like `int` or `float`.
     */
    template <typename INDEX, typename SCALAR>
    static inline SEXP spmat_to_s4(const spmat<INDEX, SCALAR> &s)
    {
      SEXP s4_class, s4;
      SEXP s4_i, s4_p, s4_Dim, s4_Dimnames, s4_x, s4_factors;
      
      const INDEX m = s.nrows();
      const INDEX n = s.ncols();
      const INDEX nnz = s.get_nnz();
      const INDEX p_len = n+1;
      
      PROTECT(s4_i = allocVector(INTSXP, nnz));
      arraytools::copy(nnz, s.index_ptr(), INTEGER(s4_i));
      
      PROTECT(s4_p = allocVector(INTSXP, p_len));
      arraytools::copy(p_len, s.col_ptr(), INTEGER(s4_p));
      
      PROTECT(s4_Dim = allocVector(INTSXP, 2));
      INTEGER(s4_Dim)[0] = (int) m;
      INTEGER(s4_Dim)[1] = (int) n;
      
      PROTECT(s4_Dimnames = allocVector(VECSXP, 2));
      SET_VECTOR_ELT(s4_Dimnames, 0, R_NilValue);
      SET_VECTOR_ELT(s4_Dimnames, 1, R_NilValue);
      
      PROTECT(s4_x = allocVector(REALSXP, nnz));
      arraytools::copy(nnz, s.data_ptr(), REAL(s4_x));
      
      PROTECT(s4_factors = allocVector(VECSXP, 0));
      
      PROTECT(s4_class = MAKE_CLASS("dgCMatrix"));
      PROTECT(s4 = NEW_OBJECT(s4_class));
      SET_SLOT(s4, install("i"), s4_i);
      SET_SLOT(s4, install("p"), s4_p);
      SET_SLOT(s4, install("Dim"), s4_Dim);
      SET_SLOT(s4, install("Dimnames"), s4_Dimnames);
      SET_SLOT(s4, install("x"), s4_x);
      SET_SLOT(s4, install("factors"), s4_factors);
      
      UNPROTECT(8);
      return s4;
    }
  }
  
  
  
  namespace internal
  {
    namespace get
    {
      template <typename INDEX, typename SCALAR>
      static inline void dim(const SEXP x, INDEX *m, INDEX *n)
      {
        spar::internal::sexp::get_dim_from_s4(x, m, n);
      }
      
      
      
      template <typename INDEX, typename SCALAR>
      static inline void col(const INDEX j, const SEXP x, spvec<INDEX, SCALAR> &s)
      {
        spar::conv::s4col_to_spvec(j, x, s);
      }
      
      
      
      template <typename INDEX, typename SCALAR>
      static inline INDEX max_col_nnz(const SEXP x)
      {
        INDEX m, n;
        spar::internal::sexp::get_dim_from_s4(x, &m, &n);
        
        SEXP P = spar::internal::sexp::get_p_from_s4(x);
        
        INDEX max_nnz = 0;
        for (INDEX col=0; col<n+1; col++)
        {
          INDEX col_nnz = spar::internal::sexp::get_col_len_from_s4(col, P);
          if (col_nnz > max_nnz)
            max_nnz = col_nnz;
        }
        
        return max_nnz;
      }
    }
  }
}


#endif
