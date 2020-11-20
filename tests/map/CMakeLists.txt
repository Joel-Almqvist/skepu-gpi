
add_executable(starpu_map starpu_map.cpp)
target_link_libraries(starpu_map
	PRIVATE
		catch2_main SkePU::SkePU MPI::MPI_CXX PkgConfig::STARPU OpenMP::OpenMP_CXX)

add_test(NAME starpu_map_seq
	COMMAND starpu_map)
add_mpi_test(starpu_map_par starpu_map)

add_executable(starpu_mapreduce starpu_mapreduce.cpp)
target_link_libraries(starpu_mapreduce
	PRIVATE
		catch2_main SkePU::SkePU MPI::MPI_CXX PkgConfig::STARPU OpenMP::OpenMP_CXX)

add_test(NAME starpu_mapreduce_seq
	COMMAND starpu_mapreduce)
add_mpi_test(starpu_mapreduce_par starpu_mapreduce)
skepu_add_executable(mvmult_cpu_test
	SKEPUSRC mvmult.cpp)
target_link_libraries(mvmult_cpu_test
	PRIVATE catch2_main)
add_test(mvmult_cpu mvmult_cpu_test)

skepu_add_executable(mvmult_openmp_test
	OpenMP
	SKEPUSRC mvmult.cpp)
target_link_libraries(mvmult_openmp_test
	PRIVATE catch2_main)
add_test(mvmult_openmp mvmult_openmp_test)

skepu_add_executable(mvmult_cuda_test
	CUDA
	SKEPUSRC mvmult.cpp)
target_link_libraries(mvmult_cuda_test
	PRIVATE catch2_main)
add_test(mvmult_cuda mvmult_cuda_test)

skepu_add_executable(mvmult_opencl_test
	OpenCL
	SKEPUSRC mvmult.cpp)
target_link_libraries(mvmult_opencl_test
	PRIVATE catch2_main)
add_test(mvmult_opencl mvmult_opencl_test)
