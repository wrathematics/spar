#include <string>

#include <gen/gen.hpp>
#include <spar.hpp>
#include <reduce.hpp>

#include "timer.hpp"

#define EARLY_EXIT -1
#define BAD_FLAG 1

using INDEX = uint32_t;
using SCALAR = int;
using MAT = spmat<INDEX, SCALAR>;


typedef struct opts_t
{ 
  bool print_header;
  bool approx;
  bool allreduce;
  INDEX n;
  float prop_dense;
} opts_t;

static inline int process_flags(int rank, int argc, char **argv, opts_t *opts)
{
  char c;
  opts->print_header = false;
  opts->approx = false;
  opts->allreduce = 0;
  opts->n = 5000;
  opts->prop_dense = 0.001;
  
  while ((c = getopt(argc, argv, "dar:n:p:h")) != -1)
  {
    if (c == 'd')
      opts->print_header = true;
    else if (c == 'a')
      opts->approx = true;
    else if (c == 'r')
      opts->allreduce = atoi(optarg);
    else if (c == 'n')
      opts->n = atoi(optarg);
    else if (c == 'p')
      opts->prop_dense = atof(optarg);
    else if (c == 'h')
    {
      if (rank == 0)
      {
        printf("Usage\n");
        printf("  mpirun -np $NRANKS ./reduce_rand\n");
        printf("Options:\n");
        printf("  -r\tTODO...\n");
      }
      
      return EARLY_EXIT;
    }
    else if (c == '?')
      return BAD_FLAG;
  }
  
  return 0;
}



int main(int argc, char **argv)
{
  opts_t opts;
  timer t;
  
  spar::mpi::init();
  int rank = spar::mpi::get_rank();
  
  int check = process_flags(rank, argc, argv, &opts);
  if (check == EARLY_EXIT || check == BAD_FLAG)
  {
    spar::mpi::finalize();
    return check;
  }
  
  
  const uint32_t seed = 1234 + rank;
  
  int root = opts.allreduce ? spar::mpi::REDUCE_TO_ALL : 0;
  INDEX n = opts.n;
  float prop_dense = opts.prop_dense;
  
  if (opts.print_header && rank == 0)
  {
    printf("seed,root,n,prop_dense,bytes_index,bytes_scalar,");
    printf("nnz_local,len_local,time_gen,");
    printf("nnz,len,time_reduce\n");
  }
  
  if (rank == 0)
    printf("%d,%d,%d,%f,%d,%d,", seed, root, n, prop_dense, (int)sizeof(INDEX), (int)sizeof(SCALAR));
  
  t.start();
  auto x = spar::gen::rand<INDEX, SCALAR>(seed, prop_dense, n, n, !opts.approx);
  t.stop();
  spar::mpi::barrier();
  if (rank == 0)
    printf("%d,%d,%f,", x.get_nnz(), x.get_len(), t.elapsed());
  
  t.start(true);
  auto y = spar::reduce::gather<MAT, INDEX, SCALAR>(root, x);
  t.stop();
  if (rank == 0)
    printf("%d,%d,%f\n", y.get_nnz(), y.get_len(), t.elapsed());
  
  spar::mpi::finalize();
  return 0;
}
