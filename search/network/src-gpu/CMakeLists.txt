#########################
#                       #
#         CMake         #
#                       #
#########################

# The minimum version of CMake necessary to build this project
cmake_minimum_required (VERSION 3.1.0)

# The name of our project
project (gwsearch-gpu)

# Adding source code files according to configuration
set (HDRS)
set (SRCS)
set (KRNS)

list (APPEND HDRS ${PROJECT_SOURCE_DIR}/auxi.h)
list (APPEND HDRS ${PROJECT_SOURCE_DIR}/floats.h)
list (APPEND HDRS ${PROJECT_SOURCE_DIR}/init.h)
list (APPEND HDRS ${PROJECT_SOURCE_DIR}/jobcore.h)
list (APPEND HDRS ${PROJECT_SOURCE_DIR}/settings.h)
list (APPEND HDRS ${PROJECT_SOURCE_DIR}/spline_z.h)
list (APPEND HDRS ${PROJECT_SOURCE_DIR}/struct.h)
list (APPEND HDRS ${PROJECT_SOURCE_DIR}/timer.h)
list (APPEND HDRS ${PROJECT_SOURCE_DIR}/CL/util.h)

list (APPEND SRCS ${PROJECT_SOURCE_DIR}/auxi.c)
list (APPEND SRCS ${PROJECT_SOURCE_DIR}/main.c)
list (APPEND SRCS ${PROJECT_SOURCE_DIR}/settings.c)
list (APPEND SRCS ${PROJECT_SOURCE_DIR}/timer.c)
list (APPEND SRCS ${PROJECT_SOURCE_DIR}/init.c)
list (APPEND SRCS ${PROJECT_SOURCE_DIR}/jobcore.c)
list (APPEND SRCS ${PROJECT_SOURCE_DIR}/spline_z.c)
list (APPEND SRCS ${PROJECT_SOURCE_DIR}/CL/util.c)

list (APPEND KRNS ${PROJECT_SOURCE_DIR}/floats.hcl)
list (APPEND KRNS ${PROJECT_SOURCE_DIR}/kernels.hcl)
list (APPEND KRNS ${PROJECT_SOURCE_DIR}/kernels.cl)

# Variable to hold ALL files to build and be visible in IDE
set (BUILD ${HDRS} ${SRCS} ${KRNS})

# Create filters for IDEs
source_group ("Headers" FILES ${HDRS})
source_group ("Sources" FILES ${SRCS})
source_group ("Kernels" FILES ${KRNS})

# Specify executable sources
add_executable (${PROJECT_NAME} ${BUILD})

# Append our project's include directory to the "#include <>" paths
target_include_directories (${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories (${PROJECT_NAME} PUBLIC ${OpenCL_INCLUDE_DIRS})
target_include_directories (${PROJECT_NAME} PUBLIC ${CLFFT_INCLUDE_DIRS})
target_include_directories (${PROJECT_NAME} PUBLIC ${CLBLAS_INCLUDE_DIRS})

# Link dependant libraries
target_link_libraries (${PROJECT_NAME} clFFT)
target_link_libraries (${PROJECT_NAME} clBLAS)

# Create project groups for IDEs
set_target_properties (${PROJECT_NAME} PROPERTIES FOLDER "GPU")

#
# Set compiler-specific variables
#
  # If compiling with Microsoft Visual C
  if (CMAKE_C_COMPILER_ID MATCHES "MSVC")

    # Generate static library of missing posix functions
    add_library (getopt SHARED ${PROJECT_SOURCE_DIR}/posix/getopt.c)
    target_include_directories (getopt PUBLIC ${PROJECT_SOURCE_DIR})
    
    # Link it to our executable
    target_link_libraries (${PROJECT_NAME} getopt)
    
    # Set definition to export symbols into library
    target_compile_definitions (getopt PRIVATE EXPORTS_GETOPT)
    
    target_compile_definitions (${PROJECT_NAME} PRIVATE _CRT_SECURE_NO_WARNINGS)
    
    # Create project groups for IDEs
    set_target_properties (${PROJECT_NAME} PROPERTIES FOLDER "GPU")

	# Set debugger options
	set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
	set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_COMMAND_ARGUEMENTS "${POLGRAW_SEARCH_NETWORK_ROOT}/../../testdata/test-data-network-injection -ident 205 -label J0322+0441 -band 000 -fpo 1391.3 -r ${POLGRAW_SEARCH_NETWORK_ROOT}/../../testdata/range-J0322+0441.txt")

  endif (CMAKE_C_COMPILER_ID MATCHES "MSVC")

#########################
#                       #
#         CTest         #
#                       #
#########################

# Add unit test
if (BUILD_TESTS)

  add_test(NAME ${PROJECT_NAME}
           COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_NAME} ${POLGRAW_SEARCH_NETWORK_ROOT}/../../testdata/test-data-network-injection -ident 205 -label J0322+0441 -band 000 -fpo 1391.3 -r ${POLGRAW_SEARCH_NETWORK_ROOT}/../../testdata/range-J0322+0441.txt")
		
endif (BUILD_TESTS)
