cmake_minimum_required(VERSION 2.8.6 FATAL_ERROR)
 
PROJECT(cloudmr)
  find_package(ITK 4.9)
  include(${ITK_USE_FILE})
  set(ITK_LIBRARIES ${ITK_LIBRARIES})



# Output directories.

if(WIN32)
  # Some libraries (e.g. vxl libs) have no dllexport markup, so we can
  # build full shared libraries only with the GNU toolchain.
  if(CMAKE_COMPILER_IS_GNUCXX)
    # CMake adds --enable-all-exports on Cygwin (since Cygwin is
    # supposed to be UNIX-like), but we need to add it explicitly for
    # a native windows build with the MinGW tools.
    if(MINGW)
      set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS
        "-shared -Wl,--export-all-symbols -Wl,--enable-auto-import")
      set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS
        "-shared -Wl,--export-all-symbols -Wl,--enable-auto-import")
      set(CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-auto-import")
    endif(MINGW)
  else(CMAKE_COMPILER_IS_GNUCXX)
   if(BUILD_SHARED_LIBS)
     set(ITK_LIBRARY_BUILD_TYPE "SHARED")
   else(BUILD_SHARED_LIBS)
     set(ITK_LIBRARY_BUILD_TYPE "STATIC")
   endif(BUILD_SHARED_LIBS)
   set(BUILD_SHARED_LIBS OFF)
  endif(CMAKE_COMPILER_IS_GNUCXX)
endif(WIN32)

 add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-std=c++11> $<$<COMPILE_LANGUAGE:C>:-std=c99>)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

   
   
set (CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fpermissive -fopenmp")

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -fopenmp")


# can we read ismrmrd?:)
ADD_SUBDIRECTORY(Tests/ReadingTest)
# yes we can!:)


# can we reconstruct image and SNR?
ADD_SUBDIRECTORY(Reconstructors/2DRecon)
# yes we can!:) RSS and B1

# can we calulate MR and PMR?
ADD_SUBDIRECTORY(Reconstructors/2DMR)
# yes we can!:)
