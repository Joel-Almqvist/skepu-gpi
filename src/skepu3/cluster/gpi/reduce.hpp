#ifndef REDUCE_HPP
#define REDUCE_HPP

#include <matrix.hpp>
#include <type_traits>
#include <numeric>
#include <cmath>

#include <GASPI.h>


namespace skepu{

  template<typename ReduceFunc>
  class Reduce1D{
  private:
    ReduceFunc rfunc;
  public:

    Reduce1D(ReduceFunc rfunc) : rfunc{rfunc} {};

     template<typename Container>
     typename Container::value_type operator()(Container& cont){
       using T = typename Container::value_type;

       if(is_skepu_container<Container>::value){

         // TODO Find a better solution than a barrier.
         // Currently it prevents multiple operations from modifying
         // the communication and container segments and needs to be at
         // the start of all functions which use these.
         gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

         gaspi_notification_id_t first_notify_id;


         T local_sum = rfunc(((T*) cont.cont_seg_ptr)[0], ((T*) cont.cont_seg_ptr)[1]);

         for(int i = 2; i < cont.local_size; i++){
           local_sum = rfunc(local_sum, ((T*) cont.cont_seg_ptr)[i]);
         }

         ((T*) cont.comm_seg_ptr)[0] = local_sum;

         int iterations = std::ceil(std::log2(cont.nr_nodes));

         bool received = true;
         int step;

         for(int i = 0; i < iterations; i++){

           step = pow(2, i);

           // Do work only if we received a value or if it is the first iteration
           if(received || i == 0){

             if(cont.rank % (step * 2) == step - 1){
               // Send

                gaspi_write_notify(cont.comm_segment_id,
                    0,
                    cont.rank + step,
                    cont.comm_segment_id + step,
                    (i + 1) * sizeof(T),
                    sizeof(T),
                    i,
                    123,
                    cont.queue,
                    GASPI_BLOCK);

               received = false;
             }
             else if(cont.rank % (step * 2) == (step * 2) - 1){
               // Receive

                gaspi_notify_waitsome(
                  cont.comm_segment_id,
                  i,
                  1,
                  &first_notify_id,
                  GASPI_BLOCK);


                ((T*) cont.comm_seg_ptr)[0] = rfunc(((T*) cont.comm_seg_ptr)[0],
                      ((T*) cont.comm_seg_ptr)[i + 1]);
             }
             else{
               // Do nothing
             }
           }
         }


         // Distribute the reduces value
         for(int i = 0; i < iterations; i++){

           step = pow(2, i);

           if(cont.rank > (cont.nr_nodes - 1) - step){

             if(cont.rank - step >= 0){

              gaspi_write_notify(cont.comm_segment_id,
                  0,
                  cont.rank - step,
                  cont.comm_segment_id - step,
                  0,
                  sizeof(T),
                  i + iterations,
                  123,
                  cont.queue,
                  GASPI_BLOCK);

            }
           }

           else if(cont.rank > (cont.nr_nodes - 1) - 2 * step){
               // receive
               gaspi_notify_waitsome(
                 cont.comm_segment_id,
                 i + iterations,
                 1,
                 &first_notify_id,
                 GASPI_BLOCK);
           }
         }

         return ((T*) cont.comm_seg_ptr)[0];
       }
       else{
         std::cout << "ERROR Non Skepu container\n";
         return typename Container::value_type{};
       }
     }

     // Should take in a backend type
     void setBackend(){}

     // Need to be implemented
     void setReduceMode(){};
  };


  // Template deduction for classes are not allowed in c++11
  // This solves this problem
  template<typename ReduceFunc>
  Reduce1D<ReduceFunc> Reduce(ReduceFunc rfunc){
    return Reduce1D<ReduceFunc>{rfunc};
  }

} // end of namespace skepu
#endif // REDUCE_HPP
