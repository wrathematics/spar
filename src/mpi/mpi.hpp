// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_MPI_MPI_H
#define SPAR_MPI_MPI_H
#pragma once


#include "defs.hpp"
#include "err.hpp"
#include "utils.hpp"


namespace spar
{
  namespace mpi
  {
    static inline void init()
    {
      MPI_Init(NULL, NULL);
    }
    
    static inline void finalize()
    {
      MPI_Finalize();
    }
    
    static inline int get_rank(MPI_Comm comm=MPI_COMM_WORLD)
    {
      int r;
      MPI_Comm_rank(comm, &r);
      return r;
    }
    
    static inline int get_size(MPI_Comm comm=MPI_COMM_WORLD)
    {
      int s;
      MPI_Comm_size(comm, &s);
      return s;
    }
    
    
    
    template <typename T>
    void reduce(int root, void *sendbuf, T *recvbuf, int count, MPI_Op op,
      MPI_Comm comm=MPI_COMM_WORLD)
    {
      int ret;
      
      const MPI_Datatype mpi_type = utils::mpi_type_lookup(*recvbuf);
      
      if (root == REDUCE_TO_ALL)
        ret = MPI_Allreduce(sendbuf, recvbuf, count, mpi_type, op, comm);
      else
        ret = MPI_Reduce(sendbuf, recvbuf, count, mpi_type, op, root, comm);
      
      err::check_MPI_ret(ret);
    }
    
    
    
    template <typename S, typename T>
    void gatherv(int root, const S *sendbuf, int sendcount, T *recvbuf,
      const int *recvcounts, const int *displs, MPI_Comm comm=MPI_COMM_WORLD)
    {
      int ret;
      
      const MPI_Datatype mpi_type_send = utils::mpi_type_lookup(*sendbuf);
      const MPI_Datatype mpi_type_recv = utils::mpi_type_lookup(*recvbuf);
      
      if (root == REDUCE_TO_ALL)
      {
        ret = MPI_Allgatherv(sendbuf, sendcount, mpi_type_send, recvbuf,
          recvcounts, displs, mpi_type_recv, comm);
      }
      else
      {
        ret = MPI_Gatherv(sendbuf, sendcount, mpi_type_send, recvbuf,
          recvcounts, displs, mpi_type_recv, root, comm);
      }
      
      err::check_MPI_ret(ret);
    }
    
    
    
    template <typename S, typename T>
    void gather(int root, const S *sendbuf, int sendcount, T *recvbuf,
      int recvcount, MPI_Comm comm=MPI_COMM_WORLD)
    {
      int ret;
      
      const MPI_Datatype mpi_type_send = utils::mpi_type_lookup(*sendbuf);
      const MPI_Datatype mpi_type_recv = utils::mpi_type_lookup(*recvbuf);
      
      if (root == REDUCE_TO_ALL)
      {
        ret = MPI_Allgather(sendbuf, sendcount, mpi_type_send, recvbuf,
          recvcount, mpi_type_recv, comm);
      }
      else
      {
        ret = MPI_Gather(sendbuf, sendcount, mpi_type_send, recvbuf,
          recvcount, mpi_type_recv, root, comm);
      }
      
      err::check_MPI_ret(ret);
    }
  }
}


#endif
