# Main CMakeLists.txt file to build library "linalg"

cmake_minimum_required(VERSION 3.20)

# set the project name
project(linalg LANGUAGES Fortran VERSION 1.0)

#####################################################################

# set default
if (NOT DEFINED BUILD_TYPE)
  set (BUILD_TYPE "RELEASE")
endif()

set (CMAKE_Fortran_FLAGS_BASE "-fPIC -ffree-line-length-none")
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_BASE} -O5")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_BASE} -g -C -Wall -fcheck=all -O")
set (CMAKE_Fortran_FLAGS_RelWithIEEE "${CMAKE_Fortran_FLAGS_RELEASE} -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans")

# set build type (DEBUG/RELEASE)
if ( ${BUILD_TYPE} MATCHES "DEBUG" )
  message("** Compiling in debug mode")
  set ( CMAKE_BUILD_TYPE Debug )
elseif ( ${BUILD_TYPE} MATCHES "RELEASE"   )
  message("** Compiling in release mode")
  set ( CMAKE_BUILD_TYPE RELEASE )
elseif ( ${BUILD_TYPE} MATCHES "IEEE"   )
  message("** Compiling in release with IEEE mode")
  set ( CMAKE_BUILD_TYPE RelWithIEEE )
else ()
  set ( CMAKE_BUILD_TYPE user )
  set ( CMAKE_Fortran_FLAGS_user  "${CMAKE_Fortran_FLAGS_BASE} ${BUILD_TYPE}" )
  message("** Compiling in user mode:  ${CMAKE_Fortran_FLAGS_user}")
endif()

#####################################################################

# Search using LAPACK_DIR or LAPACK_ROOT
if(LAPACK_DIR OR LAPACK_ROOT)
   message(VERBOSE "Directory path for Lapack: ${LAPACK_DIR}")
   message(VERBOSE "Directory path for Lapack: ${LAPACK_ROOT}")
   find_package(LAPACK REQUIRED PATHS ${LAPACK_DIR} ${LAPACK_ROOT} NO_DEFAULT_PATH)
   if (LAPACK_FOUND)
      message(STATUS "Lapack library: ${LAPACK_LIBRARIES}")
      set(LAPACK_BLAS_LIBRARIES ${LAPACK_LIBRARIES})
   else()
      MESSAGE( SEND_ERROR  "***********************************")
      MESSAGE( SEND_ERROR  "* LAPACK/BLAS LIBRARIES NOT FOUND *")
      MESSAGE( FATAL_ERROR "***********************************")
   endif()

# Search using OPENBLAS_DIR or OPENBLAS_ROOT
elseif (OpenBLAS_DIR OR OpenBLAS_ROOT)
   message(STATUS "Directory path for OpenBLAS: ${OpenBLAS_DIR}")
   message(STATUS "Directory path for OpenBLAS: ${OpenBLAS_ROOT}")
   find_package(OpenBLAS REQUIRED PATHS ${OpenBLAS_DIR} ${OpenBLAS_ROOT} NO_DEFAULT_PATH)
   if (OpenBLAS_FOUND)
      message(STATUS "OpenBLAS library: ${OpenBLAS_LIBRARIES}")
      set(LAPACK_BLAS_LIBRARIES ${OpenBLAS_LIBRARIES})
   else()
      MESSAGE( SEND_ERROR  "***********************************")
      MESSAGE( SEND_ERROR  "* LAPACK/BLAS LIBRARIES NOT FOUND *")
      MESSAGE( FATAL_ERROR "***********************************")
   endif()

else()

   # Search lapack library using default paths
   find_package(LAPACK REQUIRED)
   if (LAPACK_FOUND)
      message(VERBOSE "Lapack libraries: ${LAPACK_LIBRARIES}")
      list(GET LAPACK_LIBRARIES 0 LAPACK_LIB)
      message(STATUS "Lapack library: ${LAPACK_LIB}")
   else()
      MESSAGE( SEND_ERROR  "***************************************************")
      MESSAGE( SEND_ERROR  "* LAPACK LIBRARIES NOT FOUND **********************")
      MESSAGE( SEND_ERROR  "* Pass the full path of the directory containg it: ")
      MESSAGE( SEND_ERROR  "* -DLAPACK_DIR=/path/to/LapackConfig.cmake")
      MESSAGE( SEND_ERROR  "* -DLAPACK_ROOT=/path/to/lapackDirectory")
      MESSAGE( FATAL_ERROR "***************************************************")
   endif()

   # Search blas library using default paths
   find_package(BLAS REQUIRED)
   if (BLAS_FOUND)
      message(VERBOSE "Blas libraries: ${BLAS_LIBRARIES}")
      list(GET BLAS_LIBRARIES 0 BLAS_LIB)
      message(STATUS "Blas library: ${BLAS_LIB}")
   else()
      MESSAGE( SEND_ERROR  "***************************************************")
      MESSAGE( SEND_ERROR  "* BLAS LIBRARIES NOT FOUND ************************")
      MESSAGE( SEND_ERROR  "* Pass the full path of the directory containg it: ")
      MESSAGE( SEND_ERROR  "* -DOpenBLAS_DIR=/path/to/LapackConfig.cmake")
      MESSAGE( SEND_ERROR  "* -DOpenBLAS_ROOT=/path/to/lapackDirectory")
      MESSAGE( FATAL_ERROR "***************************************************")
   endif()

   set(LAPACK_BLAS_LIBRARIES "${LAPACK_LIBRARIES};${BLAS_LIBRARIES}")
endif()

#####################################################################

# Search "globals" library
if (NOT TARGET "globals")

   # Try to find "globals" library using
   # -Dglobals_DIR, -Dglobals_ROOT, -DGLOBALS_DIR and -DGLOBALS_ROOT
   find_package(globals PATHS ${GLOBALS_ROOT} ${GLOBALS_DIR})

   if (NOT globals_FOUND)

      # Assume that source code for globals library
      # is in the same directory of linalg library
      set(GLOBALS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../globals)

      # Assume the code is already compiled
      if (EXISTS ${GLOBALS_SRC} AND EXISTS ${GLOBALS_SRC}/build)
         find_package(globals REQUIRED
            PATHS ${GLOBALS_SRC}/build/
            NO_DEFAULT_PATH)

      # The code is not compiled
      elseif(EXISTS ${GLOBALS_SRC} AND NOT EXISTS ${GLOBALS_SRC}/build)
         add_subdirectory(${GLOBALS_SRC} ${GLOBALS_SRC}/build EXCLUDE_FROM_ALL)

      else()
         message(FATAL_ERROR "** Library 'globals' NOT FOUND!")

      endif()

   endif()

endif()

#####################################################################

# Build external libraries
add_subdirectory(external EXCLUDE_FROM_ALL)

# Build linalg library
add_subdirectory(libsrcs)
