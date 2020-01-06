#include "catch/catch.hpp"

#include <spar.hpp>



TEMPLATE_PRODUCT_TEST_CASE("construct", "[spmat]", spmat, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  const int m = 10;
  const int n = 8;
  const int len = 5;
  TestType x(m, n, len);
  
  using INDEX = decltype(+*x.index_ptr());
  
  REQUIRE( x.nrows() == (INDEX)m );
  REQUIRE( x.ncols() == (INDEX)n );
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 0 );
}



TEMPLATE_PRODUCT_TEST_CASE("insert", "[spmat]", spmat, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  const int m = 10;
  const int n = 8;
  const int len = 5;
  TestType x(m, n, len);
  
  using INDEX = decltype(+*x.index_ptr());
  using SCALAR = decltype(+*x.data_ptr());
  
  spvec<INDEX, SCALAR> s(3);
  s.insert(3, 1);
  s.insert(1, 2);
  
  x.insert(2, s);
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 2 );
  
  auto *X = x.data_ptr();
  REQUIRE( X[0] == 2 );
  REQUIRE( X[1] == 1 );
}



TEMPLATE_PRODUCT_TEST_CASE("resize", "[spmat]", spmat, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  const int m = 10;
  const int n = 8;
  int len = 2;
  TestType x(m, n, len);
  
  using INDEX = decltype(+*x.index_ptr());
  using SCALAR = decltype(+*x.data_ptr());
  
  spvec<INDEX, SCALAR> s(3);
  s.insert(1, 1);
  x.insert(2, s);
  
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 1 );
  
  len = 5;
  x.resize(len);
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 1 );
}
