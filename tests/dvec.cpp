#include "catch/catch.hpp"

#include <dvec.hpp>



TEMPLATE_PRODUCT_TEST_CASE("construct", "[dvec]", dvec, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  const int len = 5;
  TestType x(len);
  
  using INDEX = decltype(+x.get_nnz());
  
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 0 );
}



TEMPLATE_PRODUCT_TEST_CASE("insert", "[dvec]", dvec, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  const int len = 5;
  TestType x(len);
  x.insert(3, 1);
  x.insert(1, 2);
  
  using INDEX = decltype(+x.get_nnz());
  
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 2 );
  
  auto *X = x.data_ptr();
  REQUIRE( X[1] == 2 );
  REQUIRE( X[3] == 1 );
}



TEMPLATE_PRODUCT_TEST_CASE("zero", "[dvec]", dvec, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  const int len = 5;
  TestType x(len);
  
  using INDEX = decltype(+x.get_nnz());
  
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 0 );
  
  x.insert(3, 1);
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 1 );
  
  x.zero();
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 0 );
}



TEMPLATE_PRODUCT_TEST_CASE("resize", "[dvec]", dvec, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  int len = 2;
  TestType x(len);
  
  using INDEX = decltype(+x.get_nnz());
  
  x.insert(1, 1);
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 1 );
  
  len = 5;
  x.resize(len);
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 1 );
  
  auto *X = x.data_ptr();
  REQUIRE( X[4] == 0 );
  REQUIRE( X[1] == 1 );
}
