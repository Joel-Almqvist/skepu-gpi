#ifndef MAP_HPP
#define AMP_HPP

#include <matrix.hpp>
#include <type_traits>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <deque>

#include <GASPI.h>


namespace skepu{

  template<typename Function, int nr_args>
  class Map1D{
  private:
    Function func;



    /*
    * Fetches all values within the global index range:
    * [dest_cont.start_i, dest_cont.end_i ] which are stored within "ranks"
    * part of container 1 and 2.
    *
    * The values are handled in one of the 3 following ways:
    *
    * 1 - Apply func() to c1_val and c2_val and store it in dest_cont
    * 2 - Apply func() to cx_val and its matching orphaned value
    * 3 - Store the value in dest_cont and add its index to orphans
    */
    template<typename T>
    void work_on_rank(Matrix<T>& dest_cont, Matrix<T>& cont1, Matrix<T>& cont2,
      int rank, std::deque<int>& orphans){

        // std::cout << "Rank " << dest_cont.rank << " is working with local op "
        // << dest_cont.op_nr << " and remote op " << dest_cont.vclock[rank]
        // << std::endl;


      int c1_last_elem_global_index = rank != cont1.nr_nodes - 1 ?
       (rank + 1) * cont1.step -1 :
       cont1.global_size - 1;

      int c1_first_elem_global_index = rank * cont1.step;

      // Are the reads limited by dest_cont or do we want the whole of c1?
      int c1_read_from = std::max(dest_cont.start_i, c1_first_elem_global_index);
      int c1_to = std::min(dest_cont.end_i, c1_last_elem_global_index);


      int c2_last_elem_global_index = rank != cont2.nr_nodes - 1 ?
       (rank + 1) * cont2.step - 1 :
       cont2.global_size - 1;

      int c2_first_elem_global_index = rank * cont2.step;

      int c2_read_from = std::max(dest_cont.start_i, c2_first_elem_global_index);
      int c2_to = std::min(dest_cont.end_i, c2_last_elem_global_index);


      int transfer_size = dest_cont.COMM_BUFFER_NR_ELEMS / 2;



      //std::cout << "Local Rank " << dest_cont.rank << " c1_step: " << cont1.step
      //<< "  c2_step: " << cont2.step << std::endl;




      // We may not want any values one of the two containers of this rank
      if(c1_to < dest_cont.start_i){
        c1_to = -2;
      }

      if(c2_to < dest_cont.start_i){
        c2_to = -2;
      }

      // These are correct
      /*
      std::cout << "Rank " << dest_cont.rank << " owns in dest_cont: "
      << dest_cont.start_i << ", " << dest_cont.end_i <<
      " c1: " << cont1.start_i << ", " << cont1.end_i << " c2: " << cont2.start_i
      << ", " << cont2.end_i << "  and is trying to get c1: "
      << c1_read_from << ", " << c1_to << " c2: " << c2_read_from << ", " << c2_to
      << "  from rank " << rank << std::endl;
      */

      // Start at -1 so we can increment at start to make the while more readable
      int t = -1;

      int c1_send_up_to;
      int c2_send_up_to;

      int lowest_shared_i = std::max(c1_read_from, c2_read_from);

      // Globally indexed
      int c1_rec_from = -1;
      int c1_rec_to = -2;

      // Globally indexed
      int c2_rec_from = -1;
      int c2_rec_to = -2;


      bool c1_has_unread = true;
      bool c2_has_unread = true;

      int c1_sent_elems = 0;
      int c2_sent_elems = 0;



      // std::cout << "Rank " << dest_cont.rank << " Owns " << dest_cont.start_i << ", " << dest_cont.end_i
      // << " in dest cont  ";
      // std::cout << " c1 read: " << c1_read_from << ", " << c1_to;
      // std::cout <<" c2 read: " << c2_read_from << ", " << c2_to << std::endl;


      while(c1_has_unread || c2_has_unread){

        ++t;
        c1_has_unread = c1_sent_elems <= c1_to - c1_read_from;
        c2_has_unread = c2_sent_elems <= c2_to - c2_read_from;

        if(false && t == 1){
          std::cout << "Second iteration by " << dest_cont.rank <<
          " c1_sent_elems = " << c1_sent_elems << " c2_sent_elems = " <<
          c2_sent_elems <<           std::endl;
        }


        if(c1_has_unread){
          c1_rec_from = t * transfer_size + c1_read_from;
          c1_rec_to = std::min(c1_read_from + transfer_size * (t + 1) - 1, c1_to);

          c1_sent_elems += c1_rec_to - c1_rec_from + 1;

          dest_cont.read_range(c1_rec_from, c1_rec_to, 0, cont1);
        }

        if(c2_has_unread){
          c2_rec_from = t * transfer_size + c2_read_from;
          c2_rec_to = std::min(c2_read_from + transfer_size * (t + 1) - 1, c2_to);

          c2_sent_elems += c2_rec_to - c2_rec_from + 1;


          dest_cont.read_range(c2_rec_from, c2_rec_to, transfer_size * sizeof(T), cont2);
        }

        if(!c1_has_unread && !c2_has_unread){
          //std::cout << "Stop at iteration " << t << " by rank " << dest_cont.rank
          //<< std::endl;
          // Early break, we won't be transfering any more values
          break;
        }


        bool c1_has_val;
        bool c2_has_val;

        bool c1_is_orphan;
        bool c2_is_orphan;

        bool matching_values;

        std::deque<int>::iterator c1_it;
        std::deque<int>::iterator c2_it;

        T* store_at = (T*) dest_cont.cont_seg_ptr;
        T* comm_ptr = (T*) dest_cont.comm_seg_ptr;


        // Handle all values within the communication buffer
        for(int i = 0; i < transfer_size; i++){

          //std::cout << "Orphans size: " << orphans.size() << std::endl;

          matching_values = false;
          c1_has_val = c1_rec_from + i <= c1_rec_to;
          c2_has_val = c2_rec_from + i <= c2_rec_to;

          if(!c1_has_val && !c2_has_val){
            //std::cout << "Stop at iteration " << i << " by rank " << dest_cont.rank
            //<< std::endl;
            // Break early if we have no more values
            break;
          }


          if(c1_rec_from == c2_rec_from && c1_has_val && c2_has_val){
            //std::cout << "Matching value at i = " << i << std::endl;

            // If both values are on the same index we can save some operations
            // by not orphaning the values.
            matching_values = true;
            store_at[c1_rec_from - dest_cont.start_i + i]
            = func(comm_ptr[i], comm_ptr[i + transfer_size]);
          }

          if(!matching_values && c1_has_val){

            c1_it = std::find(orphans.begin(), orphans.end(), c1_rec_from + i);
            c1_is_orphan = c1_it == orphans.end();


            if(c1_is_orphan){
              /*
              std::cout << "received orphan for c1 at i = " << i
              << " with val " << comm_ptr[i] << std::endl;
              */
              orphans.push_back(c1_rec_from + i);
              store_at[c1_rec_from - dest_cont.start_i + i] = comm_ptr[i];


            }

            else{
              store_at[c1_rec_from - dest_cont.start_i + i] =
                func(
                  store_at[c1_rec_from - dest_cont.start_i + i],
                  comm_ptr[i]
                );

                /*
                std::cout << "received matched value for c1 at i = " << i << std::endl;
                std::cout << "Matching "  << store_at[c1_rec_from - dest_cont.start_i + i]
                << " with " << comm_ptr[i] << std::endl;
                */

              orphans.erase(c1_it);
            }
          }



          if(!matching_values && c2_has_val){

            c2_it = std::find(orphans.begin(), orphans.end(), c2_rec_from + i);
            c2_is_orphan = c2_it == orphans.end();

            if(c2_is_orphan){
              //std::cout << "received orphan for c2 at i = " << i << std::endl;
              orphans.push_back(c2_rec_from + i);
              store_at[c2_rec_from - dest_cont.start_i + i] = comm_ptr[i + transfer_size];
            }

            else{
              /*
              std::cout << "received matched value for c2 at i = " << i << std::endl;
              std::cout << "Matching "  << store_at[c2_rec_from - dest_cont.start_i + i]
              << " with " << comm_ptr[i+transfer_size] << std::endl;
              */

              store_at[c2_rec_from - dest_cont.start_i + i] =
                func(
                  store_at[c2_rec_from - dest_cont.start_i + i],
                  comm_ptr[i+transfer_size]
                );

              orphans.erase(c2_it);
            }
          }

        } // end of for (transfer_size)

      } // end of while()
    } // end of do_work()



        // TODO add indeces for to better see which elements we are reading from
        // each container. After reading them check if there are orphaned
        // values we should combine them with.

        // 1 - Add more indeces
        // 2 - Do we already have the corresponding value? do func()
        // 3 - After this while loop, store information about orphans

        // DO WORK






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


     /*
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

          int lowest_rank_i_depend_on = std::min(cont1.get_owner(dest_cont.start_i),
            cont2.get_owner(dest_cont.start_i));

          int highest_rank_i_depend_on = std::max(cont1.get_owner(dest_cont.end_i),
            cont2.get_owner(dest_cont.end_i));

          int lowest_rank_depending_on_me = std::min(dest_cont.get_owner(cont1.start_i),
            dest_cont.get_owner(cont2.start_i));

          int highest_rank_depending_on_me = std::max(dest_cont.get_owner(cont1.end_i),
            dest_cont.get_owner(cont2.end_i));


          dest_cont.wait_ranks.clear();


          int lowest_rank_to_wait_for = std::min(lowest_rank_i_depend_on, lowest_rank_depending_on_me);
          int highest_rank_to_wait_for = std::max(highest_rank_i_depend_on, highest_rank_depending_on_me);

          for(int i = lowest_rank_to_wait_for; i <= highest_rank_to_wait_for; i++){
            dest_cont.wait_ranks.push_back(i);
          }

          dest_cont.wait_for_vclocks(dest_cont.op_nr);

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

          dest_cont.vclock[dest_cont.rank] = ++dest_cont.op_nr;

          cont1.vclock[cont1.rank] = ++cont1.op_nr;

          cont2.vclock[cont2.rank] = ++cont2.op_nr;

     } // end of () operator
     */

     template<typename DestCont, typename Cont1, typename Cont2>
      auto func_test(DestCont& dest_cont, Cont1& cont1, Cont2& cont2) ->
      decltype(
        std::declval<typename DestCont::is_skepu_container>(),
        std::declval<typename Cont1::is_skepu_container>(),
        std::declval<typename Cont2::is_skepu_container>(),

        // Check that the lambda takes two arguments
        std::declval<Function>()(std::declval<typename DestCont::value_type>(),
        std::declval<typename DestCont::value_type>()),
        std::declval<void>()) {


          using T = typename DestCont::value_type;


          // TODO what if the rank ranges are not contiguous?
          int lowest_rank_i_depend_on = std::min(cont1.get_owner(dest_cont.start_i),
            cont2.get_owner(dest_cont.start_i));

          int highest_rank_i_depend_on = std::max(cont1.get_owner(dest_cont.end_i),
            cont2.get_owner(dest_cont.end_i));



          int lowest_local_i = std::max({cont1.start_i, cont2.start_i, dest_cont.start_i});
          int highest_local_i = std::min({cont1.end_i, cont2.end_i, dest_cont.end_i});

          // Do the local work
          for(int i = lowest_local_i; i <= highest_local_i; i++){
            ((T*) dest_cont.cont_seg_ptr)[i - dest_cont.start_i] =
            func(((T*) cont1.cont_seg_ptr)[i - cont1.start_i],
              ((T*) cont2.cont_seg_ptr)[i - cont2.start_i]);
          }


          // Manually update the vclock once since there can't be any ranks
          // which are ready at this point.
          dest_cont.get_vclock(lowest_rank_i_depend_on, dest_cont.segment_id - dest_cont.rank
           + lowest_rank_i_depend_on);

         bool got_ranks_part [highest_rank_i_depend_on - lowest_rank_i_depend_on + 1] = {};
         bool work_remaining = true;

         // Orphans are unique for each destination rank, but shared for multiple
         // calls to the same rank.
         std::deque<int> orphans{};
         int j;
         while(work_remaining){

           work_remaining = false;
           j = -1;

           for(int i = lowest_rank_i_depend_on; i <= highest_rank_i_depend_on; i++){
             j++;

             if(got_ranks_part[j] == true){
               continue;
             }

             else if(dest_cont.vclock[i] < dest_cont.op_nr){

               dest_cont.get_vclock(i, dest_cont.segment_id - dest_cont.rank + i);

               if(dest_cont.vclock[i] < dest_cont.op_nr){
                 // The updated vclock is still not ready
                 work_remaining = true;
                 continue;
               }
             }

             // std::cout << "Working with bool: " <<
             // (dest_cont.vclock[i] < dest_cont.op_nr)
             // << " vclock: " << dest_cont.vclock[i] << " op_nr: "
             // <<dest_cont.op_nr << std::endl;

             // We may or may not have updated the vclock but there is work
             // to be done at this point
             got_ranks_part[j] = true;
             work_on_rank(dest_cont, cont1, cont2, i, orphans);
           } // End of for

           assert(orphans.size() == 0);

           std::this_thread::sleep_for(std::chrono::milliseconds(5));

         }




         int lowest_rank_depending_on_me = std::min(dest_cont.get_owner(cont1.start_i),
           dest_cont.get_owner(cont2.start_i));

         int highest_rank_depending_on_me = std::max(dest_cont.get_owner(cont1.end_i),
           dest_cont.get_owner(cont2.end_i));

         dest_cont.wait_ranks.clear();
         dest_cont.vclock[dest_cont.rank] = ++dest_cont.op_nr;

         /*
         std::cout << "Rank " << dest_cont.rank << " is waiting for "
         << lowest_rank_depending_on_me << ", " << highest_rank_depending_on_me
         << " with op_nr " << dest_cont.op_nr <<std::endl;
         */

         for(int i = lowest_rank_depending_on_me; i <= highest_rank_depending_on_me; i++){
           if (i != dest_cont.rank)
             dest_cont.wait_ranks.push_back(i);
         }

         dest_cont.wait_for_vclocks(dest_cont.op_nr);

         dest_cont.vclock[dest_cont.rank] = ++dest_cont.op_nr;


         // Update vclock
         //cont1.vclock[cont1.rank] = ++cont1.op_nr;
         //cont2.vclock[cont2.rank] = ++cont2.op_nr;

         /*
         std::cout << "Rank " << dest_cont.rank << " new vlock: ";
         for(int i = 0; i < dest_cont.nr_nodes; i++){
           std::cout << dest_cont.vclock[i] << ", ";
         }
         std::cout << std::endl;
         */


         // TODO Here we must wait on ranks which depend on us





     } // end of func_test



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
