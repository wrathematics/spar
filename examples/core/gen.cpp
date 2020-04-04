#include <spar.hpp>
#include <gen/gen.hpp>


int main()
{
  printf("## A banded matrix\n");
  auto x = spar::gen::banded<int, int>(2, 10, 5);
  
  x.info();
  x.print();
  x.print(true);
  
  printf("## A band-ish matrix\n");
  auto y = spar::gen::bandish<int, int>(1234, 10, 5);
  
  y.info();
  y.print();
  y.print(true);
  
  return 0;
}
