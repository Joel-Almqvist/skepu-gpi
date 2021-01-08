/*
 * This file contains code from map which has been scrapped but might be
 * repurposed later on.
 *
 * WARNING Do NOT include this file it does not compile
 *
*/



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

    // We may not want any values from one of the two containers of this rank
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

      // TODO Parallelize this:
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
