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
  namespace internal
  {
    template <class SPMAT, typename INDEX, typename SCALAR>
    static inline INDEX get_initial_len(const SPMAT &x)
    {
      INDEX len = std::max(
        (INDEX) 32,
        (INDEX) (spar::internal::get::max_col_nnz<INDEX, SCALAR>(x) * spar::internal::defs::MEM_FUDGE_ELT_FAC)
      );
      
      return len;
    }
  }
  
  /// @brief Reducers
  namespace reduce
  {
    /**
      @brief Computes a sparse matrix (all)reduce column-by-column, where each
      column is summed via a dense vector (all)reduce.
      
      @param[in] root The number of the receiving process in the case of a
      reduce, or `spar::mpi::REDUCE_TO_ALL` for an allreduce.
      @param[in] x A supported sparse matrix in CSC format.
      @param[in] comm MPI communicator.
      
      @return An spmat object. You can convert it to an Eigen or R sparse matrix
      using the library's included converters.
      
      @comm If the input matrix has `m` rows and `n` columns, there are `n`
      (all)reduces each of length `m`.
      
      @allocs Several temporary objects are constructed. Throughout, let `m`
      denote the number of rows and `n` the number of columns of the input
      sparse matrix.
        1. (all processes) `dvec<INDEX, SCALAR>` of length `m`.
        2. (all processes) `spvec<INDEX, SCALAR>`, with initial length equal to
        the largest number of non-zero elements across all the columns (called
        `len`).
        3. (root process) The return `spmat<INDEX, SCALAR>`, with initial length
        `n*len`.
      The internal sparse vector and the return sparse matrix will resize
      themselves as needed during the reduce process.
      
      @except If there is only one MPI rank, the function will throw a
      `runtime_error` exception. If a memory allocation fails, a `bad_alloc`
      exception will be thrown. If something goes wrong with any of the MPI
      operations, a `runtime_error` exception will be thrown.
      
      @tparam SPMAT should be of type `spmat<INDEX, SCALAR>`,
      `Eigen::SparseMatrix`, or R's `dgCMatrix`.
      @tparam INDEX should be some kind of fundamental indexing type, like `int`
      or `uint16_t`.
      @tparam SCALAR should be a fundamental numeric type like `int` or `float`.
     */
    template <class SPMAT, typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> dense(const int root, const SPMAT &x, MPI_Comm comm=MPI_COMM_WORLD)
    {
      mpi::err::check_size(comm);
      const bool receiving = (root == mpi::REDUCE_TO_ALL || root == mpi::get_rank(comm));
      
      INDEX m, n;
      internal::get::dim<INDEX, SCALAR>(x, &m, &n);
      
      // setup
      const INDEX len = spar::internal::get_initial_len<SPMAT, INDEX, SCALAR>(x);
      spvec<INDEX, SCALAR> a(len);
      dvec<INDEX, SCALAR> d(m);
      
      spmat<INDEX, SCALAR> s(m, n, 0);
      if (receiving)
        s.resize(len);
      
      
      // allreduce column-by-column
      for (INDEX j=0; j<n; j++)
      {
        internal::get::col<INDEX, SCALAR>(j, x, a);
        a.densify(d);
        
        if (receiving)
          mpi::reduce(root, MPI_IN_PLACE, d.data_ptr(), m, MPI_SUM, comm);
        else
          mpi::reduce(root, d.data_ptr(), d.data_ptr(), m, MPI_SUM, comm);
        
        if (receiving)
        {
          d.update_nnz();
          a.set(d);
          s.insert(j, a);
        }
      }
      
      return s;
    }
    
    
    
    /**
      @brief Computes a sparse matrix (all)reduce column-by-column, where each
      sum is computed locally after an index (all)gather and a scalar
      (all)gather.
      
      @param[in] root The number of the receiving process in the case of a
      reduce, or `spar::mpi::REDUCE_TO_ALL` for an allreduce.
      @param[in] x A supported sparse matrix in CSC format.
      @param[in] comm MPI communicator.
      
      @return An spmat object. You can convert it to an Eigen or R sparse matrix
      using the library's included converters.
      
      @comm If the input matrix has `n` columns, there are `n` iterations of
        1. allgather the number of non-zero elements
        2. if not all of the above numbers are zero, (all)gatherv the indices
        and values
      
      @allocs Several temporary objects are constructed:
        1. (all processes) `spvec<INDEX, SCALAR>`, with initial length equal to
        the largest number of non-zero elements across all the columns (called
        `len`).
        2. (all processes) Two `dvec<int, int>` vectors, each with as many
        elements as the number of MPI ranks (denot this value as `size`). 
        3. (root process) A `std::vector<INDEX>` and a `std::vector<SCALAR>`,
        and a `std::vector<std::pair<INDEX, SCALAR>>`. All three have initial
        length `len`.
        4. (root process) The return `spmat<INDEX, SCALAR>`, with initial length
        `n*len`.
      The internal sparse vector, the three `std::vector`'s, and the return
      sparse matrix will resize themselves as needed during the reduce process.
      
      @except If there is only one MPI rank, the function will throw a
      `runtime_error` exception. If a memory allocation fails, a `bad_alloc`
      exception will be thrown. If something goes wrong with any of the MPI
      operations, a `runtime_error` exception will be thrown.
      
      @tparam SPMAT should be of type `spmat<INDEX, SCALAR>`,
      `Eigen::SparseMatrix`, or R's `dgCMatrix`.
      @tparam INDEX should be some kind of fundamental indexing type, like `int`
      or `uint16_t`.
      @tparam SCALAR should be a fundamental numeric type like `int` or `float`.
     */
    template <class SPMAT, typename INDEX, typename SCALAR>
    static inline spmat<INDEX, SCALAR> gather(const int root, const SPMAT &x, MPI_Comm comm=MPI_COMM_WORLD)
    {
      mpi::err::check_size(comm);
      const bool receiving = (root == mpi::REDUCE_TO_ALL || root == mpi::get_rank(comm));
      
      INDEX m, n;
      internal::get::dim<INDEX, SCALAR>(x, &m, &n);
      
      // setup
      const INDEX len = spar::internal::get_initial_len<SPMAT, INDEX, SCALAR>(x);
      spvec<INDEX, SCALAR> a(len);
      spmat<INDEX, SCALAR> s(m, n, 0);
      
      int size = mpi::get_size(comm);
      dvec<int, int> counts(size);
      dvec<int, int> displs(size);
      displs[0] = 0;
      
      // we need vectors of indices and values for the Allgatherv, and a vector
      // of pairs for the sort/merge
      std::vector<INDEX> indices;
      std::vector<SCALAR> values;
      std::vector<std::pair<INDEX, SCALAR>> v;
      
      if (receiving)
      {
        s.resize(len);
        
        indices.resize(len);
        values.resize(len);
        v.resize(len);
      }
      
      
      // allreduce column-by-column
      for (INDEX j=0; j<n; j++)
      {
        internal::get::col<INDEX, SCALAR>(j, x, a);
        
        // get the displacements
        INDEX count_local = a.get_nnz();
        mpi::gather(mpi::REDUCE_TO_ALL, &count_local, 1, counts.data_ptr(), 1, comm);
        
        unsigned int count = counts.sum();
        
        if (count == 0)
          continue;
        else if (receiving && indices.capacity() < count)
        {
          indices.resize(count);
          values.resize(count);
          v.resize(count);
        }
        
        for (int i=1; i<displs.get_len(); i++)
          displs[i] = displs[i-1] + counts[i-1];
        
        // get all the indices/values
        mpi::gatherv(root, a.index_ptr(), a.get_nnz(), indices.data(), counts.data_ptr(), displs.data_ptr(), comm);
        mpi::gatherv(root, a.data_ptr(),  a.get_nnz(), values.data(),  counts.data_ptr(), displs.data_ptr(), comm);
        
        // add all the vectors
        if (receiving)
        {
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
          s.insert(j, a);
        }
      }
      
      return s;
    }
  }
}


#endif
