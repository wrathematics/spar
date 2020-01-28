#include <spar.hpp>


int main()
{
  const int len = 5;
  spvec<int, int> v(len);
  v.insert(3, 1);
  v.insert(1, 2);
  
  spmat<int, int> x(10, 3, len);
  x.insert(0, v);
  
  v.zero();
  v.insert(5, 3);
  v.insert(8, 4);
  // v.insert(2, 5);
  x.insert(1, v);
  
  v.zero();
  v.insert(7, 5);
  x.insert(2, v);
  
  x.info();
  x.print();
  x.print(true);
  
  return 0;
}
