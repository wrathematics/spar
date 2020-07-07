#include <gen/gen.hpp>
#include <spar.hpp>
#include <reduce.hpp>


int main()
{
  spar::mpi::init();
  int rank = spar::mpi::get_rank();
  
  using INDEX = int;
  using SCALAR = int;
  const uint32_t seed = 1234;
  
  auto x = spar::gen::rand<INDEX, SCALAR>(seed + rank, 10, 8);
  if (rank == 0)
    x.print();
  
  auto y = spar::reduce::gather<spar::spmat<INDEX, SCALAR>, INDEX, SCALAR>(0, x);
  if (rank == 0)
    y.print();
  
  spar::mpi::finalize();
  
  return 0;
}
