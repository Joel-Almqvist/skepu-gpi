#include <iostream>
#include <GASPI.h>
#include <vector>

#include <matrix.hpp>
#include <reduce.hpp>
#include <map.hpp>
//#include <filter.hpp>


int main(){



  skepu::Matrix<long> m2{3,3,2};
  skepu::Matrix<long> m3{4,4,3};

  std::array<int, 4> arr{};

  m2.wait_for_vclocks(0, arr);

  auto map = skepu::Map<2>([](long a, long b) long {
    return a*b;
  });

  //skepu::Map<float, 2>(2.2);

   auto reduce = skepu::Reduce([](long a, long b) long {
     return a + b;
   });


  //map(m1, m2, m3);

  return 0;
}
