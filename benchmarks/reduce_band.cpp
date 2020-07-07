#include "args.hpp"
#include "common.hpp"
#include "printing.hpp"
#include "timer.hpp"

#define BENCHMARK "reduce_bandish"


int main(int argc, char **argv)
{
  opts_t<INDEX> opts;
  timer t;
  
  // setup
  spar::mpi::init();
  int rank = spar::mpi::get_rank();
  int size = spar::mpi::get_size();
  spar::mpi::err::check_size(MPI_COMM_WORLD);
  
  int check = process_flags(rank, argc, argv, &opts);
  if (check == EARLY_EXIT || check == BAD_FLAG)
  {
    spar::mpi::finalize();
    return check;
  }
  
  int root = opts.allreduce ? spar::mpi::REDUCE_TO_ALL : 0;
  INDEX n = opts.n;
  
  print_header(rank, &opts);
  print_setup<INDEX, SCALAR>(rank, size, BENCHMARK, &opts);
  
  // generation
  MAT x;
  t.start();
  if (opts.approx)
    x = spar::gen::bandish<INDEX, SCALAR>(opts.seed, n, n);
  else
    x = spar::gen::banded<INDEX, SCALAR>(opts.band, n, n);
  t.stop();
  
  print_time(rank, x, t);
  spar::mpi::barrier();
  
  // reduce
  MAT y;
  t.start(true);
  if (opts.densevec)
    y = spar::reduce::dense<MAT, INDEX, SCALAR>(root, x);
  else
    y = spar::reduce::gather<MAT, INDEX, SCALAR>(root, x);
  t.stop();
  if (rank == 0)
    printf("%d,%d,%f\n", y.get_nnz(), y.get_len(), t.elapsed());
  
  print_time(rank, y, t);
  print_final(rank);
  
  spar::mpi::finalize();
  return 0;
}
