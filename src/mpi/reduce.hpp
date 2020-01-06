// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_MPI_REDUCE_H
#define SPAR_MPI_REDUCE_H
#pragma once


#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <stdexcept>

#include "defs.hpp"
#include "../spvec.hpp"


namespace spar
{
  namespace mpi
  {
    namespace err
    {
      static inline void check_MPI_ret(int ret)
      {
        if (ret != MPI_SUCCESS)
        {
          int slen;
          char s[MPI_MAX_ERROR_STRING];
          
          MPI_Error_string(ret, s, &slen);
          throw std::runtime_error(s);
        }
      }
    }
    
    
    
    template <typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> reduce_densevec(const int root, SEXP send_data, MPI_Comm comm)
    {
      int mpi_ret;
      int m, n;
      sparsehelpers::sexp::get_dim_from_s4(send_data, &m, &n);
      
      // setup
      const int len = sparsehelpers::sexp::get_col_len_from_s4(0, sparsehelpers::sexp::get_p_from_s4(send_data)) * sparsehelpers::constants::MEM_FUDGE_ELT_FAC;
      spvec<INDEX, SCALAR> a(len);
      spmat<INDEX, SCALAR> s(m, n, n*len);
      dvec<INDEX, SCALAR> d(m);
      
      const MPI_Datatype reduce_type = mpi_type_lookup(*d.data_ptr());
      if (reduce_type == MPI_DATATYPE_NULL)
        throw std::runtime_error("unknown reducer type");
      
      // allreduce column-by-column
      for (INDEX j=0; j<n; j++)
      {
        sparsehelpers::s4col_to_spvec(j, send_data, a);
        a.densify(d);
        
        if (root == spvec::mpi::defs::REDUCE_TO_ALL)
          mpi_ret = MPI_Allreduce(MPI_IN_PLACE, d.data_ptr(), m, reduce_type, MPI_SUM, comm);
        else
          mpi_ret = MPI_Reduce(MPI_IN_PLACE, d.data_ptr(), m, reduce_type, MPI_SUM, root, comm);
        
        err::check_MPI_ret(mpi_ret);
        
        d.update_nnz();
        a.set(d);
        
        int needed_space = s.insert(j, a);
        if (needed_space > 0)
        {
          s.resize((s.get_len() + needed_space) * sparsehelpers::constants::MEM_FUDGE_ELT_FAC);
          s.insert(j, a);
        }
      }
      
      return s;
    }
  }
}


#endif
