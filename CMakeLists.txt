# Main CMakeLists.txt file to build library "linalg"

cmake_minimum_required(VERSION 3.20)

# set the project name
project(linalg LANGUAGES Fortran VERSION 1.0)

#####################################################################

# Set default build type
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Default build configuration" FORCE)

# Set flags for each build type
set(CMAKE_Fortran_FLAGS_BASE "-fPIC -ffree-line-length-none")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_BASE} -O5")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_BASE} -g -DNDEBUG -C -Wall -fcheck=all -O")
set(CMAKE_Fortran_FLAGS_RELWITHIEEE "${CMAKE_Fortran_FLAGS_RELEASE} -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans")

# Check if build type is in the list of admissible build types
set(CMAKE_CONFIGURATION_TYPES Release Debug RelWithIEEE)

list(FIND CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} FOUND_MODE)
if(FOUND_MODE EQUAL -1)
   message(WARNING "Unknown CMAKE_BUILD_TYPE '${CMAKE_BUILD_TYPE}', using 'Release'.")
   set(CMAKE_BUILD_TYPE Release)
endif(FOUND_MODE EQUAL -1)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

#####################################################################

# Find lapack and blas libraries
include(cmake/lapack_and_blas.cmake)

# Find "globals" library
include(cmake/globals.cmake)

#####################################################################

# Build external libraries
add_subdirectory(external)

# Build linalg library
add_subdirectory(libsrcs)

#####################################################################

install(EXPORT linalg-targets
   FILE linalgTargets.cmake
   DESTINATION lib/cmake/linalg
   COMPONENT library
)

include(CMakePackageConfigHelpers)

# generate the config file that includes the exports
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/linalgConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/linalg"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )

# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/linalgConfigVersion.cmake"
  VERSION "${linalg_VERSION_MAJOR}.${linalg_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion
)

# install the generated configuration files
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/linalgConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/linalgConfigVersion.cmake
  DESTINATION lib/cmake/linalg
  COMPONENT library
  )

# generate the export targets for the build tree
# needs to be after the install(TARGETS) command
export(EXPORT linalg-targets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/linalgTargets.cmake"
)
