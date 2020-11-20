#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cassert>
#include <vector>
#include <iostream>
#include <GASPI.h>
#include <type_traits>

#include <utils.hpp>
#include <container.hpp>
#include<reduce.hpp>
#include <cmath>
// TODO remove include iostream (among more?)

namespace skepu{

  template<typename T>
  class Matrix : public skepu::_gpi::Container{


    // All skeletons must be added as friend classes here
    // and in all other containers. This is an ad hoc solution due to being a demo
    template<typename TT>
    friend class Reduce1D;
  private:
    int local_size;

    // Indiciates which indeces this partition handles
    int start_i;
    int end_i;


  public:
    using value_type = T;

    Matrix(){
      std::cout << "Empty constructor called\n";
    }

    Matrix(int rows, int cols){

      // Partition the matrix so that each rank receivs an even
      // amount of elements

      int step = (rows * cols) / nr_nodes;
      // The part which can not be evenly split
      int residual = (rows * cols) % nr_nodes;

      if(rank != nr_nodes - 1){
        start_i = rank * step;
        end_i = (rank + 1) * step - 1;
      }
      else{
        start_i = rank * step;
        end_i = (rank + 1) * step + residual - 1;

      }

      local_size = end_i - start_i + 1;


      // Allow for buffering 2 * log2(N) messages
      assert(gaspi_segment_create(
        comm_segment_id,
        gaspi_size_t{2 * (sizeof(int) * ((int) std::ceil(std::log2(nr_nodes))) + 1)},
        GASPI_GROUP_ALL,
        GASPI_BLOCK,
        GASPI_ALLOC_DEFAULT
      ) == GASPI_SUCCESS);



      assert(gaspi_segment_create(cont_segment_id,
        gaspi_size_t{sizeof(T) * local_size},
        GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_ALLOC_DEFAULT
      ) == GASPI_SUCCESS);


      gaspi_segment_ptr(cont_segment_id, &cont_seg_ptr);
      gaspi_segment_ptr(comm_segment_id, &comm_seg_ptr);

      gaspi_queue_create(&queue, GASPI_BLOCK);


      // TODO remove this testing initialization
      for(int i = 0; i < local_size; i++){

        ((T*) cont_seg_ptr)[i] = T{} + start_i + i;
      }
    };
};

  template<typename T>
  struct is_skepu_container<skepu::Matrix<T>> : std::true_type {};
}

#endif // MATRIX_HPP
