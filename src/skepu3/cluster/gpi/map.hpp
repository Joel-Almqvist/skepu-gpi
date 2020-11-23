#ifndef MAP_HPP
#define AMP_HPP

#include <matrix.hpp>
#include <type_traits>
#include <numeric>
#include <cmath>


#include <GASPI.h>


namespace skepu{

  template<typename Function>
  class Map1D{
  private:
    Function func;
  public:

    Map1D(Function func) : func{func} {};

     template<typename Container>
      auto operator()(Container& cont) -> decltype(std::declval<typename Container::value_type>(),
        std::declval<void>()) {

          if(is_skepu_container<Container>::value){
            using T = typename Container::value_type;

            // Prevents race conditions on the underlying container
            gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

            for(int i = 0; i < cont.local_size; i++){
              ((T*) cont.cont_seg_ptr)[i] = func(((T*) cont.cont_seg_ptr)[i]);
            }
          }

     } // end of () overload

     // Should take in a backend type
     void setBackend(){}

     // Need to be implemented
     void setReduceMode(){};
  };


  // Template deduction for classes are not allowed in c++11
  // This solves this problem
  template<typename Function>
  Map1D<Function> Map(Function func){
    return Map1D<Function>{func};
  }

} // end of namespace skepu
#endif // MAP_HPP
