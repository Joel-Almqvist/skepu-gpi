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
    template<typename TT, int>
    friend class Map1D;
    template<typename TT>
    friend class FilterClass;
  private:

    using is_skepu_container = decltype(true);

    int local_size;
    long global_size;

    // Global information about the container
    int last_partition_size;
    long last_partition_comm_offset;

    int norm_partition_size;
    long norm_partition_comm_offset;


    long unsigned comm_size;

    // This determines how many objects of type T the communication buffer
    // will be able to fit. Too small buffer may cause major issues with Filter.
    // WARNING must be even
    static const int NR_OBJECTS_IN_COMM_BUFFER = 40;

    // Indiciates which indeces this partition handles
    int start_i;
    int end_i;

    // TODO remove this, norm_partition_size == step
    int step;


    int get_owner(int index){
      return std::min((int) std::floor((float) index / step), nr_nodes - 1);
    }



    // Puts all elements from start to end (these are global indeces) into
    // the given GASPI segment. Many to one communication pattern
    void get_range(
      int start,
      int end,
      int dest_rank,
      int dest_seg_id,
      int offset,
      int notify_id
      ) {

      int first_elem = -1;
      int last_elem = -1;

      if(end >= end_i && start <= start_i){
        // Middle part of the range
        first_elem = start_i;
        last_elem = end_i;
      }
      else if(start >= start_i && start <= end_i){
        // The start of the range
        first_elem = start;
        last_elem = end <= end_i ? end : end_i;
      }
      else if(end <= end_i && end >= start_i){
        // The end of the range
        first_elem = start >= start_i ? start : start_i;
        last_elem = end;
      }

      if(last_elem != -1){

       gaspi_write_notify(segment_id,
         sizeof(T) * (first_elem % local_size),
         dest_rank,
         dest_seg_id,
         offset + sizeof(T) * (first_elem - start), // offset remote
         sizeof(T) * (last_elem + 1 - first_elem), // size
         notify_id, // notification id
         rank + 1, // notification value, not used atm
         queue,
         GASPI_BLOCK);

      }
    }


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

      last_partition_size = step + residual;
      norm_partition_size = step;

      local_size = end_i - start_i + 1;
      global_size = rows * cols;


      comm_offset = sizeof(T) * local_size;
      last_partition_comm_offset = sizeof(T) * last_partition_size;
      norm_partition_comm_offset = sizeof(T) * norm_partition_size;

      // For filter we need two ints to show how many objects are in the buffer
      comm_size = 2 * sizeof(int) + sizeof(T) * NR_OBJECTS_IN_COMM_BUFFER;


      // This is how much memory is needed for Map with 2 argument currently
      if(comm_size < sizeof(T) * local_size * 2){
        comm_size = local_size * 2 * sizeof(T);
      }

      // Guarantee that buffer is large enough for Reduce to work
      if(comm_size < 2 * (sizeof(T) * ((int) std::ceil(std::log2(nr_nodes))) + 1)){
        throw std::invalid_argument( "Communication buffer is too small. "
      "It needs to be atleast 2 * log2(#nodes)" );
      }




      assert(gaspi_segment_create(
        segment_id,
        gaspi_size_t{sizeof(T) * local_size + comm_size},
        GASPI_GROUP_ALL,
        GASPI_BLOCK,
        GASPI_ALLOC_DEFAULT
      ) == GASPI_SUCCESS);


      gaspi_segment_ptr(segment_id, &cont_seg_ptr);
      comm_seg_ptr = ((T*) cont_seg_ptr) + local_size;

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
          segment_id,
          0,
          dest_rank,
          segment_id + dest_rank - rank, // remote segment id
          sizeof(T) * (index - step * dest_rank), // Remote offset
          sizeof(T),
          rank, // Notification id
          queue,
          GASPI_BLOCK
        );


        gaspi_notification_id_t notify_id;
        gaspi_notification_t notify_val = 0;
        gaspi_notify_waitsome(
          segment_id,
          rank,
          1,
          &notify_id,
          GASPI_BLOCK
        );
        gaspi_notify_reset(segment_id, notify_id, &notify_val);

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
