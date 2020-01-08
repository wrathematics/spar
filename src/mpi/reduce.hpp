// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_MPI_REDUCE_H
#define SPAR_MPI_REDUCE_H
#pragma once


#include <stdexcept>

#include "../spar.hpp"
#include "mpi.hpp"


namespace spar
{
  namespace mpi
  {
    template <class SPMAT, typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> reduce_densevec(const int root, const SPMAT &x, MPI_Comm comm=MPI_COMM_WORLD)
    {
      INDEX m, n;
      spar::get::dim<INDEX, SCALAR>(x, &m, &n);
      
      // setup
      const INDEX len = spar::get::max_col_nnz<INDEX, SCALAR>(x) * spar::defs::MEM_FUDGE_ELT_FAC;
      spvec<INDEX, SCALAR> a(len);
      spmat<INDEX, SCALAR> s(m, n, n*len);
      dvec<INDEX, SCALAR> d(m);
      
      // allreduce column-by-column
      for (INDEX j=0; j<n; j++)
      {
        spar::get::col<INDEX, SCALAR>(j, x, a);
        a.densify(d);
        
        reduce(root, MPI_IN_PLACE, d.data_ptr(), m, MPI_SUM, comm);
        
        d.update_nnz();
        a.set(d);
        
        INDEX needed_space = s.insert(j, a);
        if (needed_space > 0)
        {
          s.resize((s.get_len() + needed_space) * spar::defs::MEM_FUDGE_ELT_FAC);
          s.insert(j, a);
        }
      }
      
      return s;
    }
  }
}


#endif
