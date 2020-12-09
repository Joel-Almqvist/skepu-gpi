#ifndef MAP_HPP
#define AMP_HPP

#include <matrix.hpp>
#include <type_traits>
#include <numeric>
#include <cmath>


#include <GASPI.h>


namespace skepu{

  template<typename Function, int nr_args>
  class Map1D{
  private:
    Function func;
  public:

    Map1D(Function func) : func{func} {};


    template<typename DestCont, typename Cont1, typename Cont2>
     auto operator()(DestCont& dest_cont, Cont1& cont1, Cont2& cont2) ->
     decltype(
       std::declval<typename DestCont::is_skepu_container>(),
       std::declval<typename Cont1::is_skepu_container>(),
       std::declval<typename Cont2::is_skepu_container>(),

       // Check that the lambda takes one argument
       std::declval<Function>()(std::declval<typename DestCont::value_type>()),

       std::declval<void>()) {
         using T = typename DestCont::value_type;

         std::cout << "One arg map not yet implemented\n";
       }


     template<typename DestCont, typename Cont1, typename Cont2>
      auto operator()(DestCont& dest_cont, Cont1& cont1, Cont2& cont2) ->
      decltype(
        std::declval<typename DestCont::is_skepu_container>(),
        std::declval<typename Cont1::is_skepu_container>(),
        std::declval<typename Cont2::is_skepu_container>(),

        // Check that the lambda takes two arguments
        std::declval<Function>()(std::declval<typename DestCont::value_type>(),
        std::declval<typename DestCont::value_type>()),

        std::declval<void>()) {
          using T = typename DestCont::value_type;


          std::array<int, dest_cont.nr_nodes> ranks{};

          int rank_from = std::min(cont1.get_owner(dest_cont.start_i),
            cont2.get_owner(dest_cont.start_i));

          int rank_to = std::max(cont1.get_owner(dest_cont.end_i),
            cont2.get_owner(dest_cont.end_i));

          // TODO change from std::array to a ptr stored
          // in the class
          for(int i = 0; i <= nr_ranks; i++){

            if(i >= rank_from && i <= rank_to){

            }
          }



          int smallest_seg_id = dest_cont.segment_id - dest_cont.rank;

          // Each call to get range builds one nodes buffer
          for(int i = 0; i < dest_cont.nr_nodes; i++){
            if(i == cont1.nr_nodes - 1){
              // The destination node contains the unevenly split elements
              cont1.get_range(
                dest_cont.step * i,
                dest_cont.global_size - 1,
                i, // Rank
                smallest_seg_id + i, // segment id
                dest_cont.last_partition_comm_offset, // offset
                dest_cont.rank + 1 // notify id
              );

              cont2.get_range(
                dest_cont.step * i,
                dest_cont.global_size - 1,
                i, // Rank
                smallest_seg_id + i, // segment id
                dest_cont.last_partition_comm_offset + sizeof(T) * dest_cont.last_partition_size, // offset
                dest_cont.nr_nodes + dest_cont.rank + 1 // notify id
              );

            }
            else{
              cont1.get_range(
                dest_cont.step * i,
                dest_cont.step * (i + 1) - 1,
                i, // Rank
                smallest_seg_id + i, // segment id
                dest_cont.norm_partition_comm_offset, // offset
                dest_cont.rank + 1 // notify id must not be 0
              );

              cont2.get_range(
                dest_cont.step * i,
                dest_cont.step * (i + 1) - 1,
                i, // Rank
                smallest_seg_id + i, // segment id
                dest_cont.norm_partition_comm_offset + sizeof(T) * dest_cont.step, // offset
                dest_cont.nr_nodes + dest_cont.rank + 1 // notify id
              );
            }
          }


          int nr_expected_notif_c1 = cont1.get_owner(dest_cont.end_i)
              - cont1.get_owner(dest_cont.start_i) + 1;

          int nr_expected_notif_c2 = cont2.get_owner(dest_cont.end_i)
              - cont2.get_owner(dest_cont.start_i) + 1;

          gaspi_notification_t notify_val = 0;
          gaspi_notification_id_t first_id;


          for(int i = 0; i < nr_expected_notif_c1; i++){
            gaspi_notify_waitsome(
              dest_cont.segment_id , //seg id
              cont1.get_owner(dest_cont.start_i) + 1, // notification begin
              nr_expected_notif_c1, // notification number
              &first_id,
              GASPI_BLOCK
            );

            gaspi_notify_reset(dest_cont.segment_id, first_id, &notify_val);

          }


          // We need two loops since the notification ids are not contiguous
          // between the two containers
          for(int i = 0; i < nr_expected_notif_c2; i++){
            gaspi_notify_waitsome(
              dest_cont.segment_id , //seg id
              cont2.nr_nodes + cont2.get_owner(dest_cont.start_i) + 1,
              nr_expected_notif_c2, // notification number
              &first_id,
              GASPI_BLOCK
            );

            gaspi_notify_reset(dest_cont.segment_id, first_id, &notify_val);

          }

          T* dest_cont_ptr = (T*) dest_cont.cont_seg_ptr;
          T* dest_comm_ptr = (T*) dest_cont.comm_seg_ptr;

          for(int i = 0; i < dest_cont.local_size; i++){
            dest_cont_ptr[i] = func(dest_comm_ptr[i], dest_comm_ptr[i+dest_cont.local_size]);
          }

          // TODO replace with vector clock
          gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

     } // end of () operator



     // Should take in a backend type
     void setBackend(){}

     // Need to be implemented
     void setReduceMode(){};
  };


  // Template deduction for classes are not allowed in c++11
  // This solves this problem
  template<int nr_args, typename Function>
  Map1D<Function, nr_args> Map(Function func){
    return Map1D<Function, nr_args>{func};
  }



} // end of namespace skepu
#endif // MAP_HPP
