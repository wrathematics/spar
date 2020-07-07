// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_BENCHMARKS_PRINTING_HPP
#define SPAR_BENCHMARKS_PRINTING_HPP
#pragma once


#include <cstdio>

#include "args.hpp"
#include "timer.hpp"


template <typename INDEX>
static inline void print_header(const int rank, const opts_t<INDEX> *opts)
{
  if (opts->print_header && rank == 0)
  {
    printf("benchmark,");
    printf("size,");
    printf("seed,");
    printf("densevec,");
    printf("root,");
    printf("n,");
    printf("prop_dense,");
    printf("bytes_index,");
    printf("bytes_scalar,");
    printf("nnz_local,");
    printf("len_local,");
    printf("time_gen,");
    printf("nnz,");
    printf("len,");
    printf("time_reduce\n");
  }
}



template <typename INDEX, typename SCALAR>
static inline void print_setup(const int rank, const int size, const int root, const char *benchmark, const opts_t<INDEX> *opts)
{
  if (rank == 0)
  {
    printf("%s,", benchmark);
    printf("%d", size);
    printf("%d,", opts->seed);
    printf("%d,", opts->densevec);
    printf("%d,", root);
    printf("%d,", opts->n);
    printf("%f,", opts->prop_dense);
    printf("%d,", (int)sizeof(INDEX));
    printf("%d,", (int)sizeof(SCALAR));
  }
}



template <typename INDEX, typename SCALAR>
static inline void print_time(const int rank, const spmat<INDEX, SCALAR> &x, const timer &t)
{
  if (rank == 0)
  {
    printf("%d,", x.get_nnz());
    printf("%d,", x.get_len());
    printf("%f,", t.elapsed());
  }
}



static inline void print_final(const int rank)
{
  if (rank == 0)
    printf("\b \b\n");
}


#endif
