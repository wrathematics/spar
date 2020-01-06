#include <spar.hpp>


int main()
{
  const int len = 5;
  spvec<int, int> x(len);
  x.insert(3, 1);
  x.insert(1, 2);
  
  x.print();
  x.print(true);
  
  return 0;
}
