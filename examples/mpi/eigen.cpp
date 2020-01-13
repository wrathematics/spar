#include <converters/eigen.hpp>
#include <gen/gen.hpp>
#include <spar.hpp>
#include <reduce.hpp>


int main()
{
  spar::mpi::init();
  int rank = spar::mpi::get_rank();
  
  using INDEX = int;
  using SCALAR = double;
  const uint32_t seed = 1234;
  
  auto x_spmat = spar::gen::gen<INDEX, SCALAR>(seed + rank, 10, 8);
  auto x = spar::conv::spmat_to_eigen(x_spmat);
  if (rank == 0)
    std::cout << Eigen::MatrixXd(x) << std::endl;
  
  auto y = spar::reduce::gather<spar::conv::eigen_t<INDEX, SCALAR>, INDEX, SCALAR>(0, x);
  if (rank == 0)
    y.print();
  
  spar::mpi::finalize();
  
  return 0;
}
