#include <iostream>
#include <GASPI.h>
#include <vector>

#include <matrix.hpp>
//#include <reduce.hpp>
#include <map.hpp>
//#include <filter.hpp>


int main(){


  skepu::Matrix<long> m1{3,3,0};
  skepu::Matrix<long> m2{10,10,2};
  skepu::Matrix<long> m3{40,40,3};


  auto map = skepu::Map<2>([](long a, long b) long {
    return a*b;
  });

  //skepu::Map<float, 2>(2.2);

  // auto map = skepu::Map<2>([](long a) long {
  //   return a;
  // });

  map(m1, m2, m3);

  return 0;
}
