// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_REDUCE_H
#define SPAR_REDUCE_H
#pragma once


#include <vector>

#include "spar.hpp"
#include "mpi/mpi.hpp"


namespace spar
{
  /// @brief Reducers
  namespace reduce
  {
    /**
      @brief Computes a sparse matrix (all)reduce column-by-column, where each
      column is summed via a dense vector (all)reduce.
      
      @param[in] root The number of the receiving process in the case of a
      reduce, or `mpi::defs::REDUCE_TO_ALL` for an allreduce.
      @param[in] x A supported sparse matrix in CSC format.
      @param[in] comm MPI communicator.
      
      @return An spmat object. You can convert it to an Eigen or R sparse matrix
      using the library's included converters.
      
      @comm If the input matrix has `m` rows and `n` columns, there are `n`
      (all)reduces each of length `m`.
      
      @allocs Several temporary objects are constructed:
        1. A dense vector of the same fundamental type as template parameter
        `SCALAR`, with as many elements as the number of rows of the input `x`.
        2. A sparse vector (class `spvec`) with the same indexing and value
        types as the template parameters `INDEX` and `SCALAR`, respectively,
        with as many elements as the largest number of non-zero elements
        across all the columns.
        3. The return sparse matrix, initially with as many elements as the
        number of columns of the input times the maximum length described in
        item 2 above.
      The internal sparse vector and the return sparse matrix will resize
      themselves as needed during the reduce process.
      
      @except If a memory allocation fails, a `bad_alloc` exception will be
      thrown. If something goes wrong with any of the MPI operations, a
      `runtime_error` exception will be thrown.
      
      @tparam SPMAT should be of type `spmat` (included in spar),
      `Eigen::SparseMatrix`, or R's `dgCMatrix`.
      @tparam INDEX should be some kind of fundamental indexing type, like `int`
      or `uint16_t`.
      @tparam SCALAR should be a fundamental numeric type like `int` or `float`.
     */
    template <class SPMAT, typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> dense(const int root, const SPMAT &x, MPI_Comm comm=MPI_COMM_WORLD)
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
        
        mpi::reduce(root, MPI_IN_PLACE, d.data_ptr(), m, MPI_SUM, comm);
        
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
    static inline spmat<INDEX, SCALAR> gather(const int root, const SPMAT &x, MPI_Comm comm=MPI_COMM_WORLD)
    {
      INDEX m, n;
      spar::get::dim<INDEX, SCALAR>(x, &m, &n);
      
      // setup
      const INDEX len = spar::get::max_col_nnz<INDEX, SCALAR>(x) * spar::defs::MEM_FUDGE_ELT_FAC;
      spvec<INDEX, SCALAR> a(len);
      spmat<INDEX, SCALAR> s(m, n, n*len);
      
      int size = mpi::get_size(comm);
      
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
        mpi::gather(mpi::defs::REDUCE_TO_ALL, &count_local, 1, counts.data_ptr(), 1, comm);
        
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
        
        mpi::gatherv(root, a.index_ptr(), a.get_nnz(), indices.data(), counts.data_ptr(), displs.data_ptr(), comm);
        mpi::gatherv(root, a.data_ptr(),  a.get_nnz(), values.data(),  counts.data_ptr(), displs.data_ptr(), comm);
        
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
