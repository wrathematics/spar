#include <catch.hpp>
#include <spar.hpp>


TEMPLATE_PRODUCT_TEST_CASE("add sparse-sparse", "[spvec]", spvec, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  const int len = 10;
  
  TestType x(len);
  x.insert(0, 1);
  x.insert(4, 1);
  x.insert(7, 1);
  
  TestType y(len);
  y.insert(1, 2);
  y.insert(3, 2);
  y.insert(4, 2);
  y.insert(6, 2);
  y.insert(7, 1);
  y.insert(12, 2);
  
  x.add(y);
  
  REQUIRE( x.get_len() == len );
  REQUIRE( x.get_nnz() == 7 );
}



TEMPLATE_PRODUCT_TEST_CASE("add sparse-dense", "[spvec]", spvec, (
  (int, int),      (int, uint32_t),      (int, double),
  (uint32_t, int), (uint32_t, uint32_t), (uint32_t, double)
))
{
  const int len = 10;
  
  TestType x(len);
  x.insert(0, 1);
  x.insert(4, 1);
  x.insert(7, 1);
  
  const int ylen = 13;
  decltype(+*x.data_ptr()) y[ylen] = {0, 2, 2, 0, 2, 0, 2, 2, 0, 0, 0, 0, 2};
  
  x.add(y, ylen);
  
  REQUIRE( x.get_len() == len );
  REQUIRE( x.get_nnz() == 7 );
}
