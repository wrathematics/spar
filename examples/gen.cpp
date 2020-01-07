#include <spar.hpp>
#include <gen/gen.hpp>


int main()
{
  auto x = spar::gen::gen<int, int>(1234, 10, 10);
  x.print();
  
  return 0;
}
