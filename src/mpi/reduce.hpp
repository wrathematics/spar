// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_MPI_REDUCE_H
#define SPAR_MPI_REDUCE_H
#pragma once


#include <vector>

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
    
    
    
    template <class SPMAT, typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> reduce_gather(const int root, const SPMAT &x, MPI_Comm comm=MPI_COMM_WORLD)
    {
      INDEX m, n;
      spar::get::dim<INDEX, SCALAR>(x, &m, &n);
      
      // setup
      const INDEX len = spar::get::max_col_nnz<INDEX, SCALAR>(x) * spar::defs::MEM_FUDGE_ELT_FAC;
      spvec<INDEX, SCALAR> a(len);
      spmat<INDEX, SCALAR> s(m, n, n*len);
      
      int size = get_size(comm);
      
      dvec<int, int> counts(size);
      dvec<int, int> displs(size);
      displs[0] = 0;
      
      // we need vectors of indices and values for the Allgatherv, and a vector
      // of pairs for the sort/merge
      std::vector<INDEX> indices(len);
      std::vector<SCALAR> values(len);
      
      typedef std::pair<INDEX, SCALAR> pmsv; // poor man's sparse vector
      std::vector<pmsv> v;
      
      
      // allreduce column-by-column
      for (INDEX j=0; j<n; j++)
      {
        spar::get::col<INDEX, SCALAR>(j, x, a);
        
        INDEX count_local = a.get_nnz();
        gather(defs::REDUCE_TO_ALL, &count_local, 1, counts.data_ptr(), 1, comm);
        
        unsigned int count = counts.sum();
        
        if (count == 0)
          continue;
        else if (indices.capacity() < count)
        {
          indices.resize(count);
          values.resize(count);
          v.resize(count);
        }
        
        for (int i=1; i<displs.get_len(); i++)
          displs[i] = displs[i-1] + counts[i-1];
        
        gatherv(root, a.index_ptr(), a.get_nnz(), indices.data(), counts.data_ptr(), displs.data_ptr(), comm);
        gatherv(root, a.data_ptr(),  a.get_nnz(), values.data(),  counts.data_ptr(), displs.data_ptr(), comm);
        
        // add all the vectors
        for (unsigned int i=0; i<count; i++)
          v[i] = std::make_pair(indices[i], values[i]);
        
        std::sort(v.begin(), v.begin()+count);
        
        INDEX nnz = 0;
        indices[0] = v[0].first;
        values[0] = v[0].second;
        for (unsigned int i=1; i<count; i++)
        {
          if (v[i].first == indices[nnz])
            values[nnz] += v[i].second;
          else
          {
            nnz++;
            
            indices[nnz] = v[i].first;
            values[nnz] = v[i].second;
          }
        }
        
        // put summed column into the return
        a.set(nnz+1, indices.data(), values.data());
        
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
