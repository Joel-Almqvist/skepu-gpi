#include <iostream>
#include <reduce.hpp>
#include <map.hpp>
#include <GASPI.h>


main(){

  auto reduce = skepu::Reduce([](long a, long b) int {
    return a*b;
  });


  skepu::Matrix<long> matrix2{4,4};

  skepu::Matrix<long> matrix{4,4,3};

  auto map = skepu::Map([](long a) long {
    return 2;
  });


    map(matrix);
    reduce(matrix);
    matrix.print();

  return 0;
}
