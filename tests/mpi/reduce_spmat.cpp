#include <catch.hpp>
#include <spar.hpp>
#include <mpi/reduce.hpp>
#include <gen/gen.hpp>

extern int rank;
extern int size;


TEMPLATE_PRODUCT_TEST_CASE("reduce_densevec", "[spmat]", spmat, (
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
  
  spvec<INDEX, SCALAR> s(3);
  
  s.zero();
  s.insert(0, 1);
  s.insert(5, 1);
  s.insert(9, 1);
  x.insert(0, s);
  
  s.zero();
  s.insert(1, 2);
  s.insert(3, 1);
  x.insert(2, s);
  
  s.zero();
  s.insert(2, 2);
  s.insert(4, 1);
  x.insert(6, s);
  
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 7 );
  
  auto y = spar::mpi::reduce_densevec<spmat<INDEX, SCALAR>, INDEX, SCALAR>(spar::mpi::defs::REDUCE_TO_ALL, x);
  REQUIRE( y.nrows() == m );
  REQUIRE( y.ncols() == n );
  
  y.get_col(2, s);
  REQUIRE( s.get(1) == (SCALAR)2*size );
  REQUIRE( s.get(3) == (SCALAR)1*size );
}
