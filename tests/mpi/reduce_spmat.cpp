#include <catch.hpp>
#include <spar.hpp>
#include <mpi/reduce.hpp>

extern int size;


TEMPLATE_PRODUCT_TEST_CASE("reduce_densevec", "[spmat]", spmat, (
  (int, int),      (int, uint32_t),      (int, double)//,
  // (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double),
  // (int16_t, int),  (int16_t, uint32_t),  (int16_t, double),
  // (uint16_t, int), (uint16_t, uint32_t), (uint16_t, double)
))
{
  const int m = 10;
  const int n = 8;
  const int len = 5;
  TestType x(m, n, len);
  
  using INDEX = decltype(x.get_nnz());
  using SCALAR = decltype(+*x.data_ptr());
  
  spvec<INDEX, SCALAR> s(3);
  s.insert(3, 1);
  s.insert(1, 2);
  
  x.insert(2, s);
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 2 );
  
  auto y = spar::mpi::reduce_densevec<spmat<INDEX, SCALAR>, INDEX, SCALAR>(spar::mpi::defs::REDUCE_TO_ALL, x);
  REQUIRE( y.nrows() == m );
  REQUIRE( y.ncols() == n );
  
  y.get_col(2, s);
  REQUIRE( s.get(1) == (SCALAR)2*size );
  REQUIRE( s.get(3) == (SCALAR)1*size );
}
