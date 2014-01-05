cmake_minimum_required(VERSION 2.8)
set(CMAKE_VERBOSE_MAKEFILE 0)
project(CheckCUDALibs)
find_package(CUDA)

#if(CUDA_FOUND)
#	set(TEST ${CMAKE_BINARY_DIR}/test)
#	set(TESTSRC ${CMAKE_CURRENT_SOURCE_DIR}/test_capability/capability.cu)
#	execute_process(COMMAND nvcc -lcuda ${TESTSRC} -o ${TEST})
#	execute_process(COMMAND ${TEST} OUTPUT_VARIABLE CUDA_ARCH)
#	message("Cuda architecture:  ${CUDA_ARCH}")
#	list(APPEND CUDA_NVCC_FLAGS --use_fast_math ${CUDA_ARCH})
#endif(CUDA_FOUND) 

list(APPEND CUDA_NVCC_FLAGS --use_fast_math -arch=sm_21 -g)
add_definitions(-DCUDA_DEV=0 -DTIMERS=1)

if(NOT CUDA_VERSION VERSION_LESS "3.2")
  # cusparse showed up in version 3.2
  find_cuda_helper_libs(cusparse)
  find_cuda_helper_libs(curand)
  
  set(CUDA_CUSPARSE_LIBRARIES ${CUDA_cusparse_LIBRARY})
  set(CUDA_CURAND_LIBRARIES ${CUDA_curand_LIBRARY})
endif()
set(CUDA_LIBRARIES ${CUDA_LIBRARIES}; ${CUDA_CUSPARSE_LIBRARIES})
message(STATUS "CUDA_LIBRARIES = ${CUDA_LIBRARIES}")

#set(CUDA_LIBRARIES ${CUDA_LIBRARIES} 
#	"libcufft.so" 
#	"libcublas.so" 
#	"libcusparse.so")


set(CUDA_SRCS
	JobNAllSky-common.cu 
	JobNAllSky.cu  
	auxi.cu 
#	fstat.cu
	gridr.cu
	search.cu
	settings.cu)

CUDA_ADD_EXECUTABLE(gwsearch-gpu ${CUDA_SRCS})

CUDA_ADD_CUFFT_TO_TARGET(gwsearch-gpu)
CUDA_ADD_CUBLAS_TO_TARGET(gwsearch-gpu)

TARGET_LINK_LIBRARIES(gwsearch-gpu "-lrt") 


#set(BUILD_SHARED_LIBS OFF)
#set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)

