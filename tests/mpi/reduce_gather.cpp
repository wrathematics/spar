#include <catch.hpp>
#include <spar.hpp>
#include <reduce.hpp>

extern int rank;
extern int size;

#include "gen.hpp"



TEMPLATE_PRODUCT_TEST_CASE("reduce_gather", "[spmat]", spar::spmat, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double),
  (int16_t, int),  (int16_t, uint32_t),  (int16_t, double),
  (uint16_t, int), (uint16_t, uint32_t), (uint16_t, double)
))
{
  const int m = 10;
  const int n = 8;
  const int len = 10;
  TestType x(m, n, len);
  
  using INDEX = decltype(x.get_nnz());
  using SCALAR = decltype(+*x.data_ptr());
  
  fill_sparse_mat(x);
  
  auto y = spar::reduce::gather<spar::spmat<INDEX, SCALAR>, INDEX, SCALAR>(spar::mpi::REDUCE_TO_ALL, x);
  REQUIRE( y.nrows() == m );
  REQUIRE( y.ncols() == n );
  
  spar::spvec<INDEX, SCALAR> s(3);
  y.get_col(2, s);
  REQUIRE( s.get(1) == (SCALAR)2*size );
  REQUIRE( s.get(3) == (SCALAR)1*size );
  
  y.get_col(5, s);
  REQUIRE( s.get(5) == (SCALAR) 1*(size-1) );
  
  // reduce to rank 0
  auto z = spar::reduce::gather<TestType, INDEX, SCALAR>(0, x);
  REQUIRE( z.nrows() == m );
  REQUIRE( z.ncols() == n );
  
  if (rank == 0)
  {
    z.get_col(2, s);
    REQUIRE( s.get(1) == (SCALAR)2*size );
    REQUIRE( s.get(3) == (SCALAR)1*size );
    
    z.get_col(5, s);
    REQUIRE( s.get(5) == (SCALAR) 1*(size-1) );
  }
}
