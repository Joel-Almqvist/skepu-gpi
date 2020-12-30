#include <iostream>
#include <GASPI.h>
#include <vector>
#include <utility>

#include <matrix.hpp>
#include <reduce.hpp>
#include <map.hpp>
#include <omp.h>
//#include <filter.hpp>


int main(){


  skepu::Matrix<long> m1{4,4,1};
  skepu::Matrix<long> m2{5,5,2};
  skepu::Matrix<long> m3{4,4,3};


  //skepu::Matrix<long> m2{5,5,2};
  // skepu::Matrix<long> m3{4,4,3};


  auto map1 = skepu::Map<2>([](long a) long {
    return a * 2;
  });


  auto map2 = skepu::Map<2>([](long a, long b) long {
    return a + b;
  });


  for(int i = 0; i < 16; i++){
    m1.set(i,i);
  }


  map2(m1, m2, m3);
  map1(m1, m1);


  m1.print();

  return 0;
}
