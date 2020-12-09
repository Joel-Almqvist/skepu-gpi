#ifndef CONTAINER_HPP
#define CONTAINER_HPP
#include <GASPI.h>
#include <iostream>
// TODO remove iostream

namespace skepu{

  namespace _gpi{

    int curr_containers = 0;
    int created_containers = 0;

    class Container{

      // This is only a temporary, non scalable solution
      friend class Reduce1D;
      friend class Map1D;
      friend class FilterClass;
    private:
      static int counter;
      //protected:


    protected:
      gaspi_rank_t rank;
      gaspi_rank_t nr_nodes;
      gaspi_queue_id_t queue;



      gaspi_pointer_t cont_seg_ptr;
      gaspi_pointer_t comm_seg_ptr;
      long comm_offset;

    ;
      unsigned long op_nr;

      Container(){
        if(curr_containers == 0){

          // WARNING This is not a good solution, the same program may call
          // multiple init/terminates. However it is unlikely and simply doing
          // this in main leaves the remaining objects in a bad state (crashes).
          gaspi_proc_init(GASPI_BLOCK);
        }

        gaspi_proc_rank(&rank);
        gaspi_proc_num(&nr_nodes);
        segment_id = created_containers * nr_nodes + rank;

        op_nr = 0;
        curr_containers++;
        created_containers++;
      }

    public:
      // TODO move to protected
      // Must be initialized by derived classes
      unsigned long* vclock;

      // TODO Move this to protected
      gaspi_segment_id_t segment_id;

      virtual ~Container(){
        gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
        curr_containers--;
        if(curr_containers == 0){
          // WARNING This is not a good solution, the same program may call
          // multiple init/terminates. However it is unlikely and simply doing
          // this in main leaves the remaining objects in a bad state (crashes).
          gaspi_proc_term(GASPI_BLOCK);
        }
      }
    };

  }

}


#endif //CONTAINER_HPP
