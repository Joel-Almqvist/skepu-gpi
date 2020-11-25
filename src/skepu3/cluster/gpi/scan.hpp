#ifndef SCAN_HPP
#define SCAN_HPP

#include <matrix.hpp>
#include <type_traits>
#include <numeric>
#include <cmath>
#include <vector>

#include <GASPI.h>

// NOTE This skeleton will fail if the scan function gives too many results.
// Atmost NR_OBJECTS_IN_COMM_BUFFER / 2 may be found within a single partition.
// The value on this design parameter is a tradeoff between memory allocation
// and user friendliness, feel free to change it (it is defined within Matrix).

namespace skepu{

  template<typename Func>
  class Scan1D{
  private:
    Func func;
  public:

    Scan1D(Func func) : func{func} {};

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

         // TODO check if all elems fit inside comm buffer

         ((int*) cont.comm_seg_ptr)[0] = nr_elems;

         // Comm segment may contain multiple types
         T* comm_seg_ptr = ((T*) ((int*) cont.comm_seg_ptr) + 1);



         auto it = std::begin(dest);
         while(it != dest.end()){
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
  Scan1D<Func> Scan(Func func){
    return Scan1D<Func>{func};
  }

} // end of namespace skepu
#endif // SCAN_HPP
