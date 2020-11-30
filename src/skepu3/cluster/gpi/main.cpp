#include <iostream>
#include <reduce.hpp>
#include <map.hpp>
#include <scan.hpp>
#include <GASPI.h>
#include <vector>


main(){

  auto reduce = skepu::Reduce([](long a, long b) int {
    return a*b;
  });


  skepu::Matrix<long> matrix2{4,4};

  skepu::Matrix<long> matrix{20,20,3};


  matrix.set(0, 2);
  matrix.set(1, 2);
  matrix.set(2, 2);




  auto scan = skepu::Scan([](long a) bool {
    return a == 3;
  });

  std::vector<int> v{};

  scan(v, matrix);

  std::cout << v.size() << std::endl;



  return 0;
}
