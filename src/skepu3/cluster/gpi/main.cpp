#include <iostream>
#include <reduce.hpp>
#include <map.hpp>
#include <filter.hpp>
#include <GASPI.h>
#include <vector>


template<typename T1, typename T2>
T1 foofunc(T1 t1, T2 t2){
  return t1 + t2;
};

void bugg_test(){
  gaspi_proc_init(GASPI_BLOCK);
  gaspi_rank_t rank;
  gaspi_rank_t nr_nodes;

  gaspi_proc_rank(&rank);
  gaspi_proc_num(&nr_nodes);

  gaspi_queue_id_t queue;
  gaspi_queue_create(&queue, GASPI_BLOCK);

  assert(gaspi_segment_create(
    rank,
    gaspi_size_t{200},
    GASPI_GROUP_ALL,
    GASPI_BLOCK,
    GASPI_ALLOC_DEFAULT
  ) == GASPI_SUCCESS);

  gaspi_pointer_t seg_ptr;
  gaspi_segment_ptr(rank, &seg_ptr);

  if(rank == 0){
    ((int*) seg_ptr)[0] = 2;
    ((int*) seg_ptr)[1] = 3;
  }

  gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

  if(rank == 0){

    gaspi_notify(
      0, // remote seg id
      0, //rank
      10, //notification id
      33, //notif val
      queue,
      GASPI_BLOCK
    );


    gaspi_notification_id_t first_id;
    gaspi_notify_waitsome(
      0 , //seg id
      10, // notification begin
      2, // notification number
      &first_id,
      GASPI_BLOCK
    );

    gaspi_notification_t notify_val = 0;
    gaspi_notify_reset(0, first_id, &notify_val);

    std::cout << "Notify val: " << (int) notify_val << std::endl;

    gaspi_notify_reset(0, first_id, &notify_val);

    std::cout << "Notify val: " << (int) notify_val << std::endl;





    // gaspi_write_notify(0, //loc seg id
    //   0, //loc offset
    //   0, //dest rank
    //   0, //seg id remote
    //   0, // offset remote
    //   sizeof(int), // size
    //   11, // notification id
    //   333,
    //   queue,
    //   GASPI_BLOCK);


    // gaspi_write_notify(0, //loc seg id
    //   sizeof(int), //loc offset
    //   1, //dest rank
    //   1, //seg id remote
    //   sizeof(int), // offset remote
    //   sizeof(int), // size
    //   10, // notification id
    //   111,
    //   queue,
    //   GASPI_BLOCK);

  }


  else if(false && rank == 1){

    gaspi_notification_id_t first_id;



    gaspi_notify_waitsome(
      0 , //seg id
      11, // notification begin
      2, // notification number
      &first_id,
      GASPI_BLOCK
    );


    gaspi_notification_t notify_val = 0;
    gaspi_notify_reset(1, first_id, &notify_val);

    std::cout << "Notify val: " << (int) notify_val << std::endl;


    // gaspi_notify_waitsome(
    //   1 , //seg id
    //   10, // notification begin
    //   2, // notification number
    //   &first_id,
    //   GASPI_BLOCK
    // );

    gaspi_notify_reset(1, first_id, &notify_val);

    std::cout << "Notify val: " << (int) notify_val << std::endl;



    std::cout << "First id: " << (int) first_id << std::endl;
    std::cout << "Have vals: " << ((int*) seg_ptr)[0] <<", "
    << ((int*) seg_ptr)[1] << std::endl;
  }
}


int main(){

  if(false){
    bugg_test();
    return 1;
  }



  skepu::Matrix<long> m1{3,3,0};
  skepu::Matrix<long> m2{4,4,2};
  skepu::Matrix<long> m3{4,4,3};

  //gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

  /*
  auto map = skepu::Map<2>([](long a, long b) int {
    return a*b;
  });
  */

  //skepu::Map<float, 2>(2.2);

  auto map = skepu::Map<2>([](long a, long b) long {
    return a*b;
  });

  gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
  map(m1, m2, m3);

  return 0;
}
