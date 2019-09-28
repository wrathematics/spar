#include "catch/catch.hpp"

#include <spvec.hpp>



TEMPLATE_PRODUCT_TEST_CASE("construct", "[spvec]", spvec, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  const int len = 5;
  TestType x(len);
  
  REQUIRE( x.get_len() == len );
  REQUIRE( x.get_nnz() == 0 );
}



TEMPLATE_PRODUCT_TEST_CASE("insert", "[spvec]", spvec, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  const int len = 5;
  TestType x(len);
  x.insert(3, 1);
  x.insert(1, 2);
  
  REQUIRE( x.get_len() == len );
  REQUIRE( x.get_nnz() == 2 );
  
  auto *I = x.index_ptr();
  REQUIRE( I[0] == 1 );
  REQUIRE( I[1] == 3 );
  
  auto *X = x.data_ptr();
  REQUIRE( X[0] == 2 );
  REQUIRE( X[1] == 1 );
}



TEMPLATE_PRODUCT_TEST_CASE("zero", "[spvec]", spvec, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  const int len = 5;
  TestType x(len);
  
  REQUIRE( x.get_len() == len );
  REQUIRE( x.get_nnz() == 0 );
  
  x.insert(3, 1);
  REQUIRE( x.get_len() == len );
  REQUIRE( x.get_nnz() == 1 );
  
  x.zero();
  REQUIRE( x.get_len() == len );
  REQUIRE( x.get_nnz() == 0 );
}



TEMPLATE_PRODUCT_TEST_CASE("resize", "[spvec]", spvec, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  int len = 2;
  TestType x(len);
  
  x.insert(3, 1);
  REQUIRE( x.get_len() == len );
  REQUIRE( x.get_nnz() == 1 );
  
  len = 5;
  x.resize(len);
  REQUIRE( x.get_len() == len );
  REQUIRE( x.get_nnz() == 1 );
  
  auto *I = x.index_ptr();
  REQUIRE( I[0] == 3 );
  REQUIRE( I[1] == 0 );
  
  auto *X = x.data_ptr();
  REQUIRE( X[0] == 1 );
  REQUIRE( X[1] == 0 );
}
