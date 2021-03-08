# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()


# This project is designed to be built outside the Insight source tree.
project(MROptimum)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})




find_package(ISMRMRD REQUIRED)
include_directories(${ISMRMRD_INCLUDE_DIR})


add_executable(2DMR main.cxx )

target_link_libraries(2DMR ${ITK_LIBRARIES} ${Boost_LIBRARIES} ${ISMRMRD_LIBRARIES})        
