#include "catch/catch.hpp"

#include <spvec.hpp>



TEMPLATE_PRODUCT_TEST_CASE("construct", "[spvec]", spvec, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  const int len = 5;
  TestType x(len);
  
  using INDEX = decltype(+*x.index_ptr());
  
  REQUIRE( x.get_len() == (INDEX)len );
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
  
  using INDEX = decltype(+*x.index_ptr());
  
  REQUIRE( x.get_len() == (INDEX)len );
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
  
  using INDEX = decltype(+*x.index_ptr());
  
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 0 );
  
  x.insert(3, 1);
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 1 );
  
  x.zero();
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 0 );
}



TEMPLATE_PRODUCT_TEST_CASE("resize", "[spvec]", spvec, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  int len = 2;
  TestType x(len);
  
  using INDEX = decltype(+*x.index_ptr());
  
  x.insert(3, 1);
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 1 );
  
  len = 5;
  x.resize(len);
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 1 );
  
  auto *I = x.index_ptr();
  REQUIRE( I[0] == 3 );
  REQUIRE( I[1] == 0 );
  
  auto *X = x.data_ptr();
  REQUIRE( X[0] == 1 );
  REQUIRE( X[1] == 0 );
}



TEMPLATE_PRODUCT_TEST_CASE("set", "[spvec]", spvec, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  const int len = 5;
  TestType x(len);
  
  using INDEX = decltype(+*x.index_ptr());
  
  x.insert(3, 1);
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == 1 );
  
  const int setlen = 2;
  decltype(+*x.index_ptr()) I_set[setlen] = {0, 4};
  decltype(+*x.data_ptr()) X_set[setlen] = {1, 1};
  x.set(setlen, I_set, X_set);
  
  REQUIRE( x.get_len() == (INDEX)len );
  REQUIRE( x.get_nnz() == (INDEX)setlen );
  
  auto *I = x.index_ptr();
  for (INDEX i=0; i<setlen; i++)
    REQUIRE( I[i] == I_set[i] );
  for (INDEX i=setlen; i<len; i++)
    REQUIRE( I[i] == 0 );
  
  auto *X = x.data_ptr();
  for (INDEX i=0; i<setlen; i++)
    REQUIRE( X[i] == X_set[i] );
  for (INDEX i=setlen; i<len; i++)
    REQUIRE( X[i] == 0 );
}
