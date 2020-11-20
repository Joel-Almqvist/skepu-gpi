#include <iostream>
#include <reduce.hpp>
#include <GASPI.h>

class init_obj{
public:
  init_obj(){
    gaspi_proc_init(GASPI_BLOCK);
  }

  ~init_obj(){
    gaspi_proc_term(GASPI_BLOCK);
  }
};

main(){
  //init_obj obj{};
  auto reduce = skepu::Reduce([](int a, int b) int {
    //std::cout << "lambda called\n";
    return a*b;
  });



  //skepu::Matrix<int>* matrix = new skepu::Matrix<int>{5,5};
  skepu::Matrix<int> matrix{5,5};


  skepu::Matrix<int> matrix2{3,4};

  skepu::Matrix<bool> matrix3{3,3};

  // TODO set matrix values
  // TODO create map
  // TODO allow for uneven partitioning of nodes
  // TODO better idea than blocking calls?

  gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

  //skepu::Matrix<int>::value_type foo{};
  int res = reduce(matrix);
  int res2 = reduce(matrix2);
  std::cout << res << std::endl;

  std::cout << res2 << std::endl;

  // bool remote_ready = false;
  //
  // std::thread t1(
  //   [](bool& flag) -> void{
  //     flag = true;
  //     //gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
  //     //ready_flag = true;
  //     },
  //   std::ref(remote_ready));
  //
  // t1.join();



  //reduce(matrix);

  //skepu::Matrix<int>* ptr = new skepu::Matrix<int>{3,3};




  /*


  * We have access to our rank
  * We have access to the number of processes

  -> When we create a list simply partition it among number of processes

  -> For a matrix though? A smarter scheme is needed, but at compilation we do
  not know how many ranks we have.





  */


  /*
  gaspi_proc_init(GASPI_BLOCK);
  gaspi_rank_t rank;
  gaspi_rank_t num;
  gaspi_proc_rank(&rank);
  gaspi_proc_num(&num);
  printf("Hello world from rank %d of %d\n",rank, num);
  gaspi_proc_term(GASPI_BLOCK);
  */

  return 0;
}
