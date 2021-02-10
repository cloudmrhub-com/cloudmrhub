# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()


# This project is designed to be built outside the Insight source tree.
project(MROptimum)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})



find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


find_package(ISMRMRD REQUIRED)
include_directories(${ISMRMRD_INCLUDE_DIR})

ADD_DEFINITIONS(-lgsl -lgslcblas -lm)



SET(LIBS ${LIBS} ${GSL_LIBRARIES})

add_executable(2DRecon main.cxx )

target_link_libraries(2DRecon ${ITK_LIBRARIES} ${Boost_LIBRARIES} ${ISMRMRD_LIBRARIES} gsl gslcblas)  
