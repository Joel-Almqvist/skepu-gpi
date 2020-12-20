#include <iostream>
#include <GASPI.h>
#include <vector>
#include <utility>

#include <matrix.hpp>
#include <reduce.hpp>
#include <map.hpp>
//#include <filter.hpp>


int main(){

  skepu::Matrix<long> m1{4,4,1};
  skepu::Matrix<long> m2{5,5,2};
  skepu::Matrix<long> m3{4,4,3};

  //skepu::Matrix<long> m2{5,5,2};
  // skepu::Matrix<long> m3{4,4,3};


  auto map = skepu::Map<2>([](long a, long b) long {
    return a+b;
  });


  for(int i = 0; i < 16; i++){
    m2.set(i,i);
  }



  map(m1, m2, m1);

  //map(m1, m1, m1);
  //map(m1, m1, m1);
  //map(m1, m1, m1);

  m1.print();


  return 0;
}
