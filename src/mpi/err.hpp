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
    }
  }
}


#endif
