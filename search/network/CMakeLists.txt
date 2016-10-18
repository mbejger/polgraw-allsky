###############################################################################
#
# Top level CMakeLists.txt file for the polgraw-allsky repo's search/network
# Hacker: Mate Ferenc Nagy-Egri
# Last modified: 2016.10.04.
#
###############################################################################

#########################
#                       #
#         CMake         #
#                       #
#########################

# The minimum version of 'cmake' necessary to build this project
cmake_minimum_required (VERSION 3.1.0)
	
# The project name and version, release numbers
project (polgraw-search)
	
# Behavioural options for the project
option (USE_OPENCL "Build OpenCL language test suite" ON)
option (USE_CUDA "Build CUDA language test suite" OFF)
option (BUILD_TESTS "Build unit tests" OFF)
	

#
# Set platform-independent variables
#
  # Find delpendencies		
  if (${USE_OPENCL})

    find_package(OpenCL REQUIRED)
    find_package(clFFT REQUIRED CONFIG)
    find_package(clBLAS REQUIRED CONFIG)

  endif (${USE_OPENCL})
  
  if (${USE_CUDA})

    find_package(CUDA REQUIRED)

  endif (${USE_CUDA})
  
  # Use solution folders in IDEs
  set_property (GLOBAL PROPERTY USE_FOLDERS ON)
  
  # Override predefined default output directories
  set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE})

  # Set curent path for later use
  set (POLGRAW_SEARCH_NETWORK_ROOT ${PROJECT_SOURCE_DIR})
	
#
# Set platform-specific variables
#
  # If compiling on Windows
  if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		
  # If compiling on Linux-like
  if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

		
#########################
#                       #
#         CTest         #
#                       #
#########################

# Enable testing
if (BUILD_TESTS)
	
  include(CTest)
		
endif (BUILD_TESTS)
	
#########################
#                       #
#         CMake         #
#                       #
#########################
	
#
# Process targets based on user-input
#
  add_subdirectory (src-gpu)
