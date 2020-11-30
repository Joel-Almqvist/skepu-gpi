#ifndef FLITER_HPP
#define FILTER_HPP

#include <matrix.hpp>
#include <type_traits>
#include <numeric>
#include <cmath>
#include <vector>

#include <GASPI.h>

namespace skepu{

  template<typename Func>
  class FilterClass{
  private:
    Func func;
  public:

    FilterClass(Func func) : func{func} {};

    template<typename STLContainer, typename SkepuContainer>
     auto operator()(STLContainer& dest, SkepuContainer& cont) ->
     decltype(
       std::declval<typename SkepuContainer::value_type>(),
       std::begin(std::declval<STLContainer>()),
       std::end(std::declval<STLContainer>()),
       std::declval<STLContainer>().push_back(std::declval<typename STLContainer::value_type>()),
       std::declval<void>()) {


       using T = typename SkepuContainer::value_type;


       // WARNING constexpr might not be available
       if constexpr(is_skepu_container<SkepuContainer>::value){

         const long unsigned COMM_BUFFER_OFFSET{
           sizeof(int) + sizeof(T) * cont.NR_OBJECTS_IN_COMM_BUFFER / 2};


         // This barrier prevents previous function calls from modifying the
         // communication buffer
        gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

         for(int i = 0; i < cont.local_size; i++){

           if(func(((T*) cont.cont_seg_ptr)[i])){
             dest.push_back(((T*) cont.cont_seg_ptr)[i]);
           }
         }

         int nr_elems = std::end(dest) - std::begin(dest);

         int max_elems_per_iteration = cont.NR_OBJECTS_IN_COMM_BUFFER / 2;
         int received_elems[cont.nr_nodes] = {0};
         received_elems[cont.rank] = nr_elems;


         ((int*) cont.comm_seg_ptr)[0] = nr_elems;

         // Comm segment may contain multiple types
         T* comm_seg_ptr = ((T*) ((int*) cont.comm_seg_ptr) + 1);


         auto it = std::begin(dest);
         auto it_end = nr_elems <= max_elems_per_iteration ? dest.end()
            : dest.begin() + max_elems_per_iteration;

         while(it != it_end){
           *comm_seg_ptr = *it;
           it++;
           comm_seg_ptr++;
         }

         // This barrier makes sure that we do not read from an incomplete barrier
         gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);


         int next_seg;

         bool remaning_iterations = true;
         bool has_wrapped = false;

         int dest_seg = cont.comm_segment_id;
         int dest_rank;
         int min_seg_id = cont.comm_segment_id - cont.rank;
         int max_seg_id = min_seg_id + cont.nr_nodes - 1;



         int* offset_ptr;
         int nr_received_obj;

         while(remaning_iterations){
           dest_seg++;

           if(dest_seg > max_seg_id){
             dest_seg = min_seg_id;
             has_wrapped = true;
           }


           if(dest_seg == cont.comm_segment_id && has_wrapped){
             remaning_iterations = false;
             break;
           }
           else if(dest_seg == cont.comm_segment_id){
             continue;
           }

           dest_rank = cont.rank + dest_seg - cont.comm_segment_id;


           gaspi_read_notify(
             cont.comm_segment_id,
             COMM_BUFFER_OFFSET, // Local offset
             dest_rank,
             dest_seg, // remote segment id
             0, // Remote offset
             COMM_BUFFER_OFFSET, // Read size
             dest_rank, // Notification id
             cont.queue,
             GASPI_BLOCK
           );


           gaspi_notification_id_t notify_id;
           gaspi_notify_waitsome(
             cont.comm_segment_id,
             dest_rank,
             1,
             &notify_id,
             GASPI_BLOCK
           );


           // Use pointer arithmetic to find how many objects where transmitted
           offset_ptr = (int*) cont.comm_seg_ptr + 1;
           offset_ptr = (int*) (((T*) offset_ptr) + cont.NR_OBJECTS_IN_COMM_BUFFER / 2 );
           nr_received_obj = *offset_ptr;

           T* rec_objs = (T*) (offset_ptr + 1);
           for(int i = 0; i < nr_received_obj; i++){
             dest.push_back(*rec_objs);
             rec_objs++;
           }

           received_elems[cont.rank] = nr_received_obj;

           if(nr_received_obj > max_elems_per_iteration){
           }

         } // end of while

         bool remaining_work = false;

         // Check locally if work remains, otherwise this node may exit early.
         for(int i = 0; i < cont.nr_nodes; i++){
           if(received_elems[i] > max_elems_per_iteration){
             remaining_work = true;
           }
          }

          if(!remaining_work){
            return;
          }

          remaining_work = false;

         // This part is executed when the found elements did not fit in the
         // communication buffer. It is notably slower due to synchronization
         // each iteration and should be avoided by increasing buffer size.
         while(true){
           gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

           for(int i = 0; i < cont.nr_nodes; i++){
             if(received_elems[i] > max_elems_per_iteration){
               received_elems[i] = received_elems[i] - max_elems_per_iteration;

               remaining_work = true;

              if(i == cont.rank){

                // Start where we previously ended
                it = it_end;
                it_end = received_elems[i] > max_elems_per_iteration ?
                  it + max_elems_per_iteration :
                  it + received_elems[i];

                // Refill our communication buffer with elements which did not fit
                ((int*) cont.comm_seg_ptr)[0] = received_elems[i];
                comm_seg_ptr = ((T*) ((int*) cont.comm_seg_ptr) + 1);

                while(it != it_end){
                  *comm_seg_ptr = *it;
                  it++;
                  comm_seg_ptr++;
                }
              }
             }
             else{
               received_elems[i] = 0;
             }
           }

           if(!remaining_work){
             return;
           }


           for(int i = 0; i < cont.nr_nodes; i++){

             if(received_elems[i] != 0 && cont.rank != i){

               gaspi_read_notify(
                 cont.comm_segment_id,
                 COMM_BUFFER_OFFSET, // Local offset
                 i,
                 cont.comm_segment_id - cont.rank + i, // remote segment id
                 0, // Remote offset
                 COMM_BUFFER_OFFSET, // Read size
                 i, // Notification id
                 cont.queue,
                 GASPI_BLOCK
               );


               gaspi_notification_id_t notify_id;
               gaspi_notify_waitsome(
                 cont.comm_segment_id,
                 i,
                 1,
                 &notify_id,
                 GASPI_BLOCK
               );

               // Use pointer arithmetic to find how many objects where transmitted
               offset_ptr = (int*) cont.comm_seg_ptr + 1;
               offset_ptr = (int*) (((T*) offset_ptr) + max_elems_per_iteration );

               received_elems[i] = *offset_ptr;

               T* rec_objs = (T*) (offset_ptr + 1);
               for(int i = 0; i < received_elems[i]; i++){
                 dest.push_back(*rec_objs);
                 rec_objs++;
               }
             }
           }
           remaining_work = false;
         }
       }
    }


     // Should take in a backend type
     void setBackend(){}

     // Need to be implemented
     void setReduceMode(){};
  };


  // Template deduction for classes are not allowed in c++11
  // This solves this problem
  template<typename Func>
  FilterClass<Func> Filter(Func func){
    return FilterClass<Func>{func};
  }

} // end of namespace skepu
#endif // FILTER_HPP
