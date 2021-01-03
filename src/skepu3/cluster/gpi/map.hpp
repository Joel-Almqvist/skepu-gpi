#ifndef MAP_HPP
#define AMP_HPP

#include <matrix.hpp>

#include <type_traits>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <deque>
#include <utility>
#include <mutex>

#include <omp.h>
#include <GASPI.h>


namespace skepu{

  template<typename Function, int nr_args>
  class Map1D{
  private:
    Function func;



    /*
    * Fetches the values at the global index [dest_cont.start_i, dest_cont.end_i]
    * of container 1 and container 2. Ignores any index which the specified
    * rank is not of the owner of.
    *
    * WARNING: This function assumes distinct containers and it will
    * cause errors otherwise.
    *
    * The values are handled in one of the 3 following ways:
    *
    * 1 - Apply func() to c1_val and c2_val and store it in dest_cont
    * 2 - Apply func() to cx_val and its matching orphaned value
    * 3 - Store the value in orphans
    */
    template<typename T>
    void apply_rank_unique(Matrix<T>& dest_cont, Matrix<T>& cont1, Matrix<T>& cont2,
      int rank, std::deque<std::tuple<unsigned int, unsigned int, T>>& orphans){

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


      // We may not want any values one of the two containers of this rank
      if(c1_to < dest_cont.start_i){
        c1_to = -2;
      }

      if(c2_to < dest_cont.start_i){
        c2_to = -2;
      }

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

      std::mutex vlock;

      // Loop while there are elements to transfer from c1 or c2
      while(c1_has_unread || c2_has_unread){

        ++t;
        c1_has_unread = c1_sent_elems <= c1_to - c1_read_from;
        c2_has_unread = c2_sent_elems <= c2_to - c2_read_from;

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
          // Early break, we won't be transfering any more values
          break;
        }


        T* store_at = (T*) dest_cont.cont_seg_ptr;
        T* comm_ptr = (T*) dest_cont.comm_seg_ptr;


        std::vector<std::tuple<unsigned int, unsigned int, T>> new_orphans{};
        #pragma omp_parallel parallel
        {
          for(int i = omp_get_thread_num(); i < transfer_size;
          i = i + omp_get_num_threads()){

            if(c1_rec_from + i > c1_rec_to){
              // Base case
              break;
            }

            // if the matching value exists in the buffer
            if(c1_rec_from + i >= c2_rec_from && c1_rec_from + i <= c2_rec_to){
              int pair_offset = transfer_size + c1_rec_from + i - c2_rec_from;

              store_at[c1_rec_from - dest_cont.start_i + i]
              = func(comm_ptr[i], comm_ptr[pair_offset]);
            }

            else{
              vlock.lock();
              new_orphans.push_back(std::tuple<unsigned int, unsigned int, T>
                {c1_rec_from + i, 0, comm_ptr[i]});
              vlock.unlock();
            }
          }


          for(int i = omp_get_thread_num(); i < transfer_size;
          i = i + omp_get_num_threads()){

            if(c2_rec_from + i > c2_rec_to){
              // Base case
              break;
            }

            if(c2_rec_from + i >= c1_rec_from && c2_rec_from + i <= c1_rec_to){
              // We have a matching pair but the previous loop already managed
              // this case
              continue;
            }
            else{
              vlock.lock();
              new_orphans.push_back(std::tuple<unsigned int, unsigned int, T>
                {c2_rec_from + i, 0, comm_ptr[transfer_size + i]});
              vlock.unlock();
            }
          }

          // Match the new orphans with the previous ones
          for(int i = omp_get_thread_num(); i < new_orphans.size();
          i = i + omp_get_num_threads()){
            auto it = std::find_if(orphans.begin(), orphans.end(), [&new_orphans, i]
            (typename std::tuple<unsigned int, unsigned int, T> tup){
              return std::get<0>(tup) ==std::get<0>(new_orphans[i]);
              }
            );

            if(it != orphans.end()){
              int index = std::get<0>(*it);
              store_at[index - dest_cont.start_i] = func(std::get<2>(*it),
              std::get<2>(new_orphans[i]));

              // Mark the two ex orphans for removal
              std::get<1>(new_orphans[i]) = 1;
              std::get<1>(*it) = 1;

            }
          }
        } // End of parallel region



        orphans.erase(
          std::remove_if(orphans.begin(), orphans.end(), []
          (typename std::tuple<unsigned int, unsigned int, T> tup){
            return std::get<1>(tup) == 1;
          }),
          orphans.end()
        );

        std::copy_if(
          new_orphans.begin(), new_orphans.end(), std::back_inserter(orphans),
          [](typename std::tuple<unsigned int, unsigned int, T> tup){
            return std::get<1>(tup) != 1;
          }
        );
      } // end of while()
    } // end of apply_rank_unique()


    template<typename T>
    void fill_remote(Matrix<T>& dest_cont, Matrix<T>& cont1, Matrix<T>& cont2,
      int rank, std::vector<std::tuple<int, T, T>>& remote){

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


        int transfer_size = dest_cont.COMM_BUFFER_NR_ELEMS;


        // We may not want any values one of the two containers of this rank
        if(c1_to < dest_cont.start_i){
          c1_to = -2;
        }

        if(c2_to < dest_cont.start_i){
          c2_to = -2;
        }

        int t = -1;
        int rec_to;
        int rec_from;
        T* comm_ptr = (T*) dest_cont.comm_seg_ptr;
        typename std::vector<std::tuple<int, T, T>>::iterator it;
        while(true){
          t++;
          rec_to = std::min(c1_read_from + (t + 1) * transfer_size, c1_to);
          rec_from = c1_read_from + t * transfer_size;

          if(rec_from > rec_to){
            break;
          }

          dest_cont.read_range(rec_from, rec_to, 0, cont1);

          for(int i = 0; i < transfer_size; i++){

            if(i > rec_to - rec_from ){
              // We are not using the full transfer_size
              break;
            }

            it = std::find_if(remote.begin(), remote.end(),
            [rec_from, i](std::tuple<int, T, T>& tup){
              return std::get<0>(tup) == rec_from + i;
            });

            if(it == remote.end()){
              remote.push_back(std::tuple<int, T, T>
                {rec_from + i, comm_ptr[i], T{}});
            }
            else{
              std::get<1>(*it) = comm_ptr[i];
            }
          }
        }

        t = -1;
        while(true){
          t++;
          rec_to = std::min(c2_read_from + (t + 1) * transfer_size, c2_to);
          rec_from = c2_read_from + t * transfer_size;

          if(rec_from > rec_to){
            break;
          }

          dest_cont.read_range(rec_from, rec_to, 0, cont2);

          for(int i = 0; i < transfer_size; i++){

            if(i > rec_to - rec_from ){
              // We are not using the full transfer_size
              break;
            }

            it = std::find_if(remote.begin(), remote.end(),
            [rec_from, i](std::tuple<int, T, T>& tup){
              return std::get<0>(tup) == rec_from + i;
            });

            if(it == remote.end()){
              remote.push_back(std::tuple<int, T, T>
                {rec_from + i, T{}, comm_ptr[i]});
            }
            else{
              std::get<2>(*it) = comm_ptr[i];
            }
          }
        }
      }


  public:

    Map1D(Function func) : func{func} {};


    /* Performs Map with 1 argument the following way:
    * 1 - Apply map to locally owned elements
    * 2 - Wait for all ranks which have elements we need to access remotely
    * 3 - Gather remote elements and apply Map on them
    * 4 - Wait until all remote ranks are done reading from us
    */
    template<typename DestCont, typename From>
     auto operator()(DestCont& dest_cont, From& from) ->
     decltype(
       std::declval<typename DestCont::is_skepu_container>(),
       std::declval<typename From::is_skepu_container>(),

       // Check that the lambda takes one argument
       std::declval<Function>()(std::declval<typename DestCont::value_type>()),

       std::declval<void>()) {
         using T = typename DestCont::value_type;

         int lowest_local_i = std::max({from.start_i, dest_cont.start_i});
         int highest_local_i = std::min({from.end_i, dest_cont.end_i});

         T* dest_ptr = (T*) dest_cont.cont_seg_ptr;
         T* from_ptr = (T*) from.cont_seg_ptr;


         // Do the local work
         #pragma omp_parallel parallel shared(dest_cont, dest_ptr, from_ptr, \
           from, lowest_local_i, highest_local_i)
         {

           for(int i = lowest_local_i + omp_get_thread_num();
             i <= highest_local_i;
             i = i + omp_get_num_threads()){
               dest_ptr[i - dest_cont.start_i] = func(from_ptr[i - from.start_i]);
           }
         }


         int lowest_rank_i_depend_on = from.get_owner(dest_cont.start_i);
         int highest_rank_i_depend_on = from.get_owner(dest_cont.end_i);

         dest_cont.wait_ranks.clear();
         for(int i = lowest_rank_i_depend_on; i <= highest_rank_i_depend_on; i++){
           if (i != dest_cont.rank)
             dest_cont.wait_ranks.push_back(i);
         }


         dest_cont.wait_for_vclocks(dest_cont.op_nr);

         int transfer_size = dest_cont.COMM_BUFFER_NR_ELEMS;

         int range1_start = dest_cont.start_i;
         int range1_end = from.start_i;

         int range2_start = from.end_i;
         int range2_end = dest_cont.end_i;

         int start;
         int end;
         int t = -1;

         T* dest_cptr = (T*) dest_cont.comm_seg_ptr;
         if(range1_end > range1_start){
           start = range1_start;
           end = start;

           while(true){
             ++t;
             start = start + t * transfer_size;

             if(start >= range1_end){
               break;
             }

             // Minues one because read_range is inclusive
             end = std::min(end + (t + 1) * transfer_size, range1_end -1);
             dest_cont.read_range(start, end, 0, from);


             for(int i = 0; i <= end - start; i++){
               dest_ptr[start + i - dest_cont.start_i] = func(dest_cptr[i]);
             }
           }
         }


         if(range2_end > range2_start){

           start = range2_start;
           end = start;
           t = -1;

           while(true){
             ++t;
             start = start + t * transfer_size;

             if(start >= range2_end){
               break;
             }

             // Minues one because read_range is inclusive
             end = std::min(end + (t + 1) * transfer_size, range2_end -1);
             dest_cont.read_range(start, end, 0, from);

             for(int i = 0; i <= end - start; i++){
               dest_ptr[start + i - dest_cont.start_i] = func(dest_cptr[i]);
             }
           }
         }

         // Indicate that the first phase of the operation is done
         dest_cont.vclock[dest_cont.rank] = ++dest_cont.op_nr;


         // We must wait for other ranks which may want to read from us
         int lowest_rank_depending_on_me = dest_cont.get_owner(from.start_i);
         int highest_rank_depending_on_me = dest_cont.get_owner(from.end_i);

         dest_cont.wait_ranks.clear();

         for(int i = lowest_rank_depending_on_me; i <= highest_rank_depending_on_me; i++){
           if (i != dest_cont.rank)
             dest_cont.wait_ranks.push_back(i);
         }

         dest_cont.wait_for_vclocks(dest_cont.op_nr);

         dest_cont.vclock[dest_cont.rank] = ++dest_cont.op_nr;

       }


     /*
     * This function calls the correct implementation of 2 args Map
     * depending on if the given containers are unique or not.
     */
     template<typename DestCont, typename Cont1, typename Cont2>
      auto operator()(DestCont& dest_cont, Cont1& cont1, Cont2& cont2) ->
      decltype(
        std::declval<typename DestCont::is_skepu_container>(),
        std::declval<typename Cont1::is_skepu_container>(),
        std::declval<typename Cont2::is_skepu_container>(),

        // Check that the lambda takes two arguments
        std::declval<Function>()(std::declval<typename DestCont::value_type>(),
        std::declval<typename DestCont::value_type>()),
        std::declval<void>()){


          if(cont1 == cont2 || cont1 == dest_cont || cont2 == dest_cont){

            apply_on_shared_conts(dest_cont, cont1, cont2);

          }
          else{
            apply_on_unique_conts(dest_cont, cont1, cont2);
          }


        }


     template<typename DestCont, typename Cont1, typename Cont2>
      void apply_on_unique_conts(DestCont& dest_cont, Cont1& cont1,
        Cont2& cont2)
      {

          using T = typename DestCont::value_type;

          int lowest_rank_i_depend_on = std::min(cont1.get_owner(dest_cont.start_i),
            cont2.get_owner(dest_cont.start_i));

          int highest_rank_i_depend_on = std::max(cont1.get_owner(dest_cont.end_i),
            cont2.get_owner(dest_cont.end_i));



          int lowest_local_i = std::max({cont1.start_i, cont2.start_i, dest_cont.start_i});
          int highest_local_i = std::min({cont1.end_i, cont2.end_i, dest_cont.end_i});

          // Do the local work
          #pragma omp_parallel parallel shared(dest_cont, cont2, cont1, lowest_local_i, highest_local_i)
          {

            for(int i = lowest_local_i + omp_get_thread_num();
              i <= highest_local_i;
              i = i + omp_get_num_threads()){
              ((T*) dest_cont.cont_seg_ptr)[i - dest_cont.start_i] =
              func(((T*) cont1.cont_seg_ptr)[i - cont1.start_i],
              ((T*) cont2.cont_seg_ptr)[i - cont2.start_i]);
            }
          }


          // Manually update the vclock once since there can't be any ranks
          // which are ready at this point.
          dest_cont.get_vclock(lowest_rank_i_depend_on, dest_cont.segment_id - dest_cont.rank
           + lowest_rank_i_depend_on);

         bool got_ranks_part [highest_rank_i_depend_on - lowest_rank_i_depend_on + 1] = {};
         bool work_remaining = true;

         // Orphans are unique for each destination rank, but shared for multiple
         // calls to the same rank.
         std::deque<std::tuple<unsigned int, unsigned int, T>> orphans{};
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

             // We may or may not have updated the vclock but there is work
             // to be done at this point
             got_ranks_part[j] = true;
             apply_rank_unique(dest_cont, cont1, cont2, i, orphans);
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


         for(int i = lowest_rank_depending_on_me; i <= highest_rank_depending_on_me; i++){
           if (i != dest_cont.rank)
             dest_cont.wait_ranks.push_back(i);
         }

         dest_cont.wait_for_vclocks(dest_cont.op_nr);

         dest_cont.vclock[dest_cont.rank] = ++dest_cont.op_nr;


     } // end of apply_on_unique_conts



     /*
     * In this case we need to have our remote data read before we work on it locally
     */
     template<typename DestCont, typename Cont1, typename Cont2>
      void apply_on_shared_conts(DestCont& dest_cont, Cont1& cont1,
        Cont2& cont2)
      {

          using T = typename DestCont::value_type;

          int lowest_rank_i_depend_on = std::min(cont1.get_owner(dest_cont.start_i),
            cont2.get_owner(dest_cont.start_i));

          int highest_rank_i_depend_on = std::max(cont1.get_owner(dest_cont.end_i),
            cont2.get_owner(dest_cont.end_i));


          int lowest_local_i = std::max({cont1.start_i, cont2.start_i, dest_cont.start_i});
          int highest_local_i = std::min({cont1.end_i, cont2.end_i, dest_cont.end_i});


         bool got_ranks_part [highest_rank_i_depend_on - lowest_rank_i_depend_on + 1] = {};
         bool work_remaining = true;

         // Contains all remote data we will read.Theoretically it may grow very
         // large, but in practice this is unlikely.
         std::vector<std::tuple<int, T, T>> remote{};


         int j;
         while(work_remaining){

           work_remaining = false;
           j = -1;

           for(int i = lowest_rank_i_depend_on; i <= highest_rank_i_depend_on; i++){
             j++;

             if(got_ranks_part[j] == true || i == dest_cont.rank){
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

             // We may or may not have updated the vclock but there is work
             // to be done at this point
             got_ranks_part[j] = true;
             fill_remote(dest_cont, cont1, cont2, i, remote);
           } // End of for

           std::this_thread::sleep_for(std::chrono::milliseconds(5));

         }


         int lowest_rank_depending_on_me = std::min(dest_cont.get_owner(cont1.start_i),
           dest_cont.get_owner(cont2.start_i));

         int highest_rank_depending_on_me = std::max(dest_cont.get_owner(cont1.end_i),
           dest_cont.get_owner(cont2.end_i));

         dest_cont.wait_ranks.clear();
         dest_cont.vclock[dest_cont.rank] = ++dest_cont.op_nr;


         for(int i = lowest_rank_depending_on_me; i <= highest_rank_depending_on_me; i++){
           if (i != dest_cont.rank)
             dest_cont.wait_ranks.push_back(i);
         }

         dest_cont.wait_for_vclocks(dest_cont.op_nr);


         T val1;
         T val2;
         int index;

         T* dest_ptr = (T*) dest_cont.cont_seg_ptr;
         T* c1_ptr = (T*) cont1.cont_seg_ptr;
         T* c2_ptr = (T*) cont2.cont_seg_ptr;

         for(int i = 0; i < remote.size(); i++){
           index = std::get<0>(remote[i]);

           if(index >= cont1.start_i && index <= cont1.end_i){
             val1 = c1_ptr[index - cont1.start_i];
             val2 = std::get<2>(remote[i]);
           }
           else if(index >= cont2.start_i && index <= cont2.end_i){
             val1 = std::get<1>(remote[i]);
             val2 = c2_ptr[index - cont2.start_i];
           }
           else{
             val1 = std::get<1>(remote[i]);
             val2 = std::get<2>(remote[i]);
           }
           dest_ptr[index - dest_cont.start_i] = func(val1, val2);
         }

         // Do the pure local work
         #pragma omp_parallel parallel shared(dest_cont, cont2, cont1, lowest_local_i, highest_local_i)
         {
           for(int i = lowest_local_i + omp_get_thread_num();
             i <= highest_local_i;
             i = i + omp_get_num_threads()){
               dest_ptr[i - dest_cont.start_i] =
               func(((T*) cont1.cont_seg_ptr)[i - cont1.start_i],
               ((T*) cont2.cont_seg_ptr)[i - cont2.start_i]);
           }
         }


         dest_cont.vclock[dest_cont.rank] = ++dest_cont.op_nr;


     } // end of apply on non unique




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
