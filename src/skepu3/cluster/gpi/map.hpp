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
        std::declval<void>()) {
          using T = typename DestCont::value_type;


          if(!std::is_same<T, long>::value){
            std::cout << "NOT SAME\n";
            return;
          }



          int nr_notifications = 1 + cont1.get_owner(dest_cont.end_i)
              - cont1.get_owner(dest_cont.start_i);


          int smallest_seg_id = dest_cont.comm_segment_id - cont1.rank;


          // Each call to get range builds one nodes buffer
          for(int i = 0; i < cont1.nr_nodes; i++){
            if(i == cont1.nr_nodes - 1){
              cont1.get_range(
                dest_cont.step * i,
                dest_cont.global_size - 1,
                i, // Rank
                smallest_seg_id + i, // segment id
                0, // offset
                dest_cont.rank + 1 // notify id
              );

              cont2.get_range(
                dest_cont.step * i,
                dest_cont.global_size - 1,
                i, // Rank
                smallest_seg_id + i, // segment id
                sizeof(T) * dest_cont.last_partition_size, // offset
                dest_cont.nr_nodes + dest_cont.rank + 1 // notify id
              );

            }
            else{
              cont1.get_range(
                dest_cont.step * i,
                dest_cont.step * (i + 1) - 1,
                i, // Rank
                smallest_seg_id + i, // segment id
                0, // offset
                dest_cont.rank + 1 // notify id must not be 0
              );

              cont2.get_range(
                dest_cont.step * i,
                dest_cont.step * (i + 1) - 1,
                i, // Rank
                smallest_seg_id + i, // segment id
                sizeof(T) * dest_cont.step, // offset
                dest_cont.nr_nodes + dest_cont.rank + 1 // notify id
              );
            }
          }



          int nr_expected_notif_c1 = cont1.get_owner(dest_cont.end_i)
              - cont1.get_owner(dest_cont.start_i) + 1;
          int nr_expected_notif_c2 = cont2.get_owner(dest_cont.end_i)
              - cont2.get_owner(dest_cont.start_i) + 1;

          int nr_notif = 0;
          gaspi_notification_t notify_val = 0;
          gaspi_notification_id_t first_id;


          for(int i = 0; i < nr_expected_notif_c1; i++){


          //   if(dest_cont.rank != 3){
          //     return;
          //   }
          //
          //   std::cout << "Owner of " << dest_cont.start_i <<
          //   " is: " << cont1.get_owner(dest_cont.start_i) << std::endl;
          //
          //
          // // TODO get_range() skickar aldrig lokala delar till buffern
          // std::cout << "Rank " << dest_cont.rank <<
          // " seg " << (int) dest_cont.comm_segment_id  << " is waiting on seg "
          // << (int) dest_cont.comm_segment_id << " for "
          // << nr_expected_notif_c1 << " notifs starting at " <<
          // cont1.get_owner(dest_cont.start_i) + 1 << std::endl;
          // return;

            gaspi_notify_waitsome(
              dest_cont.comm_segment_id , //seg id
              cont1.get_owner(dest_cont.start_i) + 1, // notification begin
              nr_expected_notif_c1, // notification number
              &first_id,
              GASPI_BLOCK
            );
            gaspi_notify_reset(dest_cont.comm_segment_id, first_id, &notify_val);

            std::cout << "Rank " << dest_cont.rank << " got notif: ("
            << (int) first_id << ", " << (int) notify_val << ")\n";
          }


          // We need two loops since the notification ids are not contiguous
          // between the two containers
          for(int i = 0; i < nr_expected_notif_c2; i++){
            gaspi_notify_waitsome(
              dest_cont.comm_segment_id , //seg id
              cont2.nr_nodes + cont2.get_owner(dest_cont.start_i) + 1,
              nr_expected_notif_c2, // notification number
              &first_id,
              GASPI_BLOCK
            );

            gaspi_notify_reset(dest_cont.comm_segment_id, first_id, &notify_val);

            std::cout << "Rank " << dest_cont.rank << " got notif: ("
            << (int) first_id << ", " << (int) notify_val << ")\n";
          }


          /*
          std::cout << "Rank " << dest_cont.rank << " waiting for " <<
          nr_notifications << " notifications, first id to wait for: "
          << cont1.get_owner(dest_cont.start_i) <<
          " has comm seg id: "<< (int) dest_cont.comm_segment_id <<std::endl;
          */

          T* dest_cont_ptr = (T*) dest_cont.cont_seg_ptr;
          T* dest_comm_ptr = (T*) dest_cont.comm_seg_ptr;
          for(int i = 0; i < dest_cont.local_size; i++){
            // std::cout << "Adding: " << dest_comm_ptr[i] << " and " <<
            // dest_comm_ptr[i+dest_cont.local_size] << std::endl;
            // std::cout << "loc size: " << dest_cont.local_size << std::endl;

            dest_cont_ptr[i] = func(dest_comm_ptr[i], dest_comm_ptr[i+dest_cont.local_size]);
          }

          // TODO
          // 1 - After we are done with this function we can NOT guarantee that
          // that the remote state is correct. Hence we can NOT read nor write
          // from remote containers.

          // Send notifications to every node after every operation, then at the
          // start of every function wait for the nodes which we will read from.

          // Weak synchronization vs hard synchronization


          /*
          At the start of all function call hold(ranks)
          which makes these rank sleep


          A - At the end of all functions call notify to all nodes which are relying
          on me.
          B - Then wait for notifies from all nodes which rely on me.

          OR replace B with a wait at the start of all functions and store the
          relevant nodes in the matrix class. More error prone but possible better
          performance

          This scheme is hindered by the fact that all containers are distributed
          equally among all nodes. Hence the difference from this and a barrier
          is negligeble in the systems current state.

          */



          // 2 - Handle small containers

          // 3 . Expand the interface to match SkePU, add iterators



          // TODO remove this debug
          for(int i = 0; i < 4; i++){
            if(dest_cont.rank == i){
              std::cout << "Rank " << i << ": ";
              for(int j = 0; j < dest_cont.local_size; j++){
                std::cout << dest_cont_ptr[j] << ", ";
              }
              std::cout << std::endl;
            }
            gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

          }


          // TODO change this
          /*
          gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
          nr_notifications = cont2.get_range(
            dest_cont.start_i,
            dest_cont.end_i,
            dest_cont.rank, // This will be the same regardless of container
            dest_cont.comm_segment_id,
            (dest_cont.end_i + 1) * sizeof(cont1.value_type)
          )
          */

          /*
          if(is_skepu_container<Container>::value){
            using T = typename Container::value_type;

            // Prevents race conditions on the underlying container
            gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

            for(int i = 0; i < cont.local_size; i++){
              ((T*) cont.cont_seg_ptr)[i] = func(((T*) cont.cont_seg_ptr)[i]);
            }
          }
          */
     } // end of () overload



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
