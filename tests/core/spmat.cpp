#include <catch.hpp>
#include <spar.hpp>


TEMPLATE_PRODUCT_TEST_CASE("construct", "[spmat]", spar::spmat, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double),
  (int16_t, int),  (int16_t, uint32_t),  (int16_t, double),
  (uint16_t, int), (uint16_t, uint32_t), (uint16_t, double)
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



TEMPLATE_PRODUCT_TEST_CASE("insert", "[spmat]", spar::spmat, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double),
  (int16_t, int),  (int16_t, uint32_t),  (int16_t, double),
  (uint16_t, int), (uint16_t, uint32_t), (uint16_t, double)
))
{
  const int m = 10;
  const int n = 8;
  const int len = 5;
  TestType x(m, n, len);
  
  using INDEX = decltype(x.get_nnz());
  using SCALAR = decltype(+*x.data_ptr());
  
  spar::spvec<INDEX, SCALAR> s(3);
  s.insert(3, 1);
  s.insert(1, 2);
  
  x.insert(2, s);
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 2 );
  
  auto *X = x.data_ptr();
  REQUIRE( X[0] == 2 );
  REQUIRE( X[1] == 1 );
}



TEMPLATE_PRODUCT_TEST_CASE("zero", "[spmat]", spar::spmat, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double),
  (int16_t, int),  (int16_t, uint32_t),  (int16_t, double),
  (uint16_t, int), (uint16_t, uint32_t), (uint16_t, double)
))
{
  const int m = 10;
  const int n = 8;
  const int len = 5;
  TestType x(m, n, len);
  
  using INDEX = decltype(x.get_nnz());
  using SCALAR = decltype(+*x.data_ptr());
  
  spar::spvec<INDEX, SCALAR> s(3);
  s.insert(3, 1);
  s.insert(1, 2);
  
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 0 );
  
  x.insert(2, s);
  
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 2 );
  
  x.zero();
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 0 );
}



TEMPLATE_PRODUCT_TEST_CASE("resize", "[spmat]", spar::spmat, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double),
  (int16_t, int),  (int16_t, uint32_t),  (int16_t, double),
  (uint16_t, int), (uint16_t, uint32_t), (uint16_t, double)
))
{
  const int m = 10;
  const int n = 8;
  int len = 2;
  TestType x(m, n, len);
  
  using INDEX = decltype(x.get_nnz());
  using SCALAR = decltype(+*x.data_ptr());
  
  spar::spvec<INDEX, SCALAR> s(3);
  s.insert(1, 1);
  x.insert(2, s);
  
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 1 );
  
  len = 5;
  x.resize(len);
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 1 );
}



TEMPLATE_PRODUCT_TEST_CASE("get_col", "[spmat]", spar::spmat, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double),
  (int16_t, int),  (int16_t, uint32_t),  (int16_t, double),
  (uint16_t, int), (uint16_t, uint32_t), (uint16_t, double)
))
{
  const int m = 10;
  const int n = 8;
  TestType x(m, n, 20);
  
  using INDEX = decltype(x.get_nnz());
  using SCALAR = decltype(+*x.data_ptr());
  
  spar::spvec<INDEX, SCALAR> s(8);
  
  s.insert(1, 1);
  s.insert(4, 2);
  x.insert(0, s);
  
  s.zero();
  s.insert(3, 1);
  s.insert(5, 2);
  s.insert(6, 3);
  x.insert(2, s);
  
  s.zero();
  s.insert(4, 1);
  s.insert(6, 2);
  s.insert(7, 3);
  s.insert(9, 4);
  x.insert(6, s);
  
  s.zero();
  s.insert(2, 1);
  x.insert(7, s);
  
  x.get_col(0, s);
  REQUIRE( s.get_nnz() == 2 );
  
  x.get_col(1, s);
  REQUIRE( s.get_nnz() == 0 );
  
  x.get_col(2, s);
  REQUIRE( s.get_nnz() == 3 );
  
  x.get_col(5, s);
  REQUIRE( s.get_nnz() == 0 );
  
  x.get_col(7, s);
  REQUIRE( s.get_nnz() == 1 );
}
