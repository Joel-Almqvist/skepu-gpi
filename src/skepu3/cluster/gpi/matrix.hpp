#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cassert>
#include <vector>
#include <iostream>
#include <GASPI.h>
#include <type_traits>

#include <utils.hpp>
#include <container.hpp>
#include <reduce.hpp>
#include <cmath>
// TODO remove include iostream (among more?)

namespace skepu{

  template<typename T>
  class Matrix : public skepu::_gpi::Container{


    // All skeletons must be added as friend classes here
    // and in all other containers. This is an ad hoc solution due to being a demo
    template<typename TT>
    friend class Reduce1D;
    template<typename TT>
    friend class Map1D;
    template<typename TT>
    friend class FilterClass;
  private:
    int local_size;
    long unsigned comm_size;

    // This determines how many objects of type T the communication buffer
    // will be able to fit. Too small buffer may cause major issues with Scan.
    // WARNING must be even
    static const int NR_OBJECTS_IN_COMM_BUFFER = 40;

    // Indiciates which indeces this partition handles
    int start_i;
    int end_i;

    int step;

  public:

    using value_type = T;

    Matrix(){
      std::cout << "Empty constructor called\n";
    }

    Matrix(int rows, int cols){
      // Partition the matrix so that each rank receivs an even
      // amount of elements
      step = (rows * cols) / nr_nodes;
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



      comm_size = 2*sizeof(int) + sizeof(T) * NR_OBJECTS_IN_COMM_BUFFER;

      // Guarantee that buffer is large enough for Reduce to work
      if(comm_size < 2 * (sizeof(T) * ((int) std::ceil(std::log2(nr_nodes))) + 1)){
        throw std::invalid_argument( "Communication buffer is too small. "
      "It needs to be atleast 2 * log2(#nodes)" );
      }


      assert(gaspi_segment_create(
        comm_segment_id,
        gaspi_size_t{comm_size},
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

    };


    Matrix(int rows, int cols, T init) : Matrix(rows, cols){
      set(init);
    }


    void set(T scalar){
      for(int i = 0; i < local_size; i ++){
        ((T*) cont_seg_ptr)[i] = scalar;
      }
    }


    // Randomly initializes an int or long matrix
    void rand(int from, int to){
      if(std::is_same<T, int>::value || std::is_same<T, long>::value){
        for(int i = 0; i < local_size; i ++){
          ((T*) cont_seg_ptr)[i] = from + std::rand() % (to - from);
        }
      }
      else{
        std::cout << "Rand only supports matrices of int or long\n";
      }
    }



    void set(int index, T value){
      if(index >= start_i && index <= end_i){
        ((T*) cont_seg_ptr)[index - start_i] = value;
      }
    }


    T get(int index){
      // Prevents unexpected changes in the communication buffer
      gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

      if(index >= start_i && index <= end_i){
        // The value is local
        return ((T*) cont_seg_ptr)[index - start_i];
      }
      else{

        // The rank holding the value
        int dest_rank = std::floor(((double) index) / step);
        if(dest_rank >= nr_nodes){
          dest_rank = nr_nodes - 1;
        }

        gaspi_read_notify(
          comm_segment_id,
          0,
          dest_rank,
          cont_segment_id + dest_rank - rank, // remote segment id
          sizeof(T) * (index - step * dest_rank), // Remote offset
          sizeof(T),
          rank, // Notification id
          queue,
          GASPI_BLOCK
        );


        gaspi_notification_id_t notify_id;
        gaspi_notify_waitsome(
          comm_segment_id,
          rank,
          1,
          &notify_id,
          GASPI_BLOCK
        );

        return ((T*) comm_seg_ptr)[0];
      }
    }



    // WARNING Only use this for debugging, it has very poor performance
    void print(){
      for(int i = 0; i < nr_nodes; i++){
        if(i == rank){
          for(int j = 0; j < local_size; j++){
            std::cout << ((T*) cont_seg_ptr)[j] << ", ";
          }
          std::cout << std::endl;
        }
        gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
      }
    }


};

  template<typename T>
  struct is_skepu_container<skepu::Matrix<T>> : std::true_type {};
}

#endif // MATRIX_HPP
