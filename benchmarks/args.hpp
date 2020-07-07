// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_BENCHMARKS_ARGS_HPP
#define SPAR_BENCHMARKS_ARGS_HPP
#pragma once


#include <unistd.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#define EARLY_EXIT -1
#define BAD_FLAG 1


template <typename INDEX>
struct opts_t
{ 
  bool print_header;
  
  bool densevec;
  bool allreduce;
  INDEX n;
  uint32_t seed;
  
  bool approx;
  float prop_dense;
  
  INDEX band;
};



template <typename INDEX>
static inline int process_flags(int rank, int argc, char **argv, opts_t<INDEX> *opts)
{
  char c;
  opts->print_header = false;
  
  opts->densevec = false;
  opts->allreduce = false;
  opts->n = 5000;
  opts->seed = 1234;
  
  opts->approx = false;
  opts->prop_dense = 0.001;
  
  opts->band = 1;
  
  while ((c = getopt(argc, argv, "davrn:p:s:b:h")) != -1)
  {
    if (c == 'd')
      opts->print_header = true;
    else if (c == 'a')
      opts->approx = true;
    else if (c == 'v')
      opts->densevec = true;
    else if (c == 'r')
      opts->allreduce = true;
    else if (c == 'n')
      opts->n = atoi(optarg);
    else if (c == 'p')
      opts->prop_dense = atof(optarg);
    else if (c == 's')
      opts->seed = atof(optarg) + rank;
    else if (c == 'b')
      opts->band = atof(optarg);
    else if (c == 'h')
    {
      if (rank == 0)
      {
        printf("Usage\n");
        printf("  mpirun -np $NRANKS ./reduce_rand\n");
        // printf("Options:\n");
        // printf("  -r\tTODO...\n");
      }
      
      return EARLY_EXIT;
    }
    else if (c == '?')
      return BAD_FLAG;
  }
  
  return 0;
}


#endif
