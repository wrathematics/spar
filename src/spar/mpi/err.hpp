// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_MPI_ERR_H
#define SPAR_MPI_ERR_H
#pragma once


#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <cstdint>
#include <typeinfo>


namespace spar
{
  namespace mpi
  {
    static inline int get_size(MPI_Comm comm);
    static inline void finalize();
    
    
    
    namespace err
    {
      static inline void check_ret(int ret)
      {
        if (ret != MPI_SUCCESS)
        {
          int slen;
          char s[MPI_MAX_ERROR_STRING];
          
          MPI_Error_string(ret, s, &slen);
          throw std::runtime_error(s);
        }
      }
      
      
      
      static inline void check_size(MPI_Comm comm)
      {
        if (spar::mpi::get_size(comm) < 2)
        {
          spar::mpi::finalize();
          throw std::runtime_error("reducer requires more than 1 MPI rank");
        }
      }
    }
  }
}


#endif
