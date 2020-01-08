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
    template <typename T>
    void reduce_inplace(int root, T *recvbuf, int count, MPI_Op op, MPI_Comm comm)
    {
      int ret;
      const MPI_Datatype mpi_type = utils::mpi_type_lookup(*recvbuf);
      
      if (root == defs::REDUCE_TO_ALL)
        ret = MPI_Allreduce(MPI_IN_PLACE, recvbuf, count, mpi_type, op, comm);
      else
        ret = MPI_Reduce(MPI_IN_PLACE, recvbuf, count, mpi_type, op, root, comm);
      
      err::check_MPI_ret(ret);
    }
    
    
    template <typename T>
    void gatherv(int root, const T *sendbuf, int sendcount, T *recvbuf,
      const int *recvcounts, const int *displs, MPI_Comm comm)
    {
      int ret;
      const MPI_Datatype mpi_type = utils::mpi_type_lookup(*recvbuf);
      
      if (root == defs::REDUCE_TO_ALL)
        ret = MPI_Allgatherv(sendbuf, sendcount, mpi_type, recvbuf, recvcounts, displs, mpi_type, comm);
      else
        ret = MPI_Gatherv(sendbuf, sendcount, mpi_type, recvbuf, recvcounts, displs, mpi_type, root, comm);
      
      err::check_MPI_ret(ret);
    }
  }
}


#endif
