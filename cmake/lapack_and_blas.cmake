# Library is specified by the user
if (LAPACK_LIBRARIES)
   if (NOT EXISTS ${LAPACK_LIBRARIES})
      message(FATAL_ERROR "Lapack library does not exists: ${LAPACK_LIBRARIES}")
   endif ()

# Find directly the library
elseif (LAPACK_LIBRARY_DIR)
   find_library(LAPACK_LIBRARIES
      NAME lapack
      PATHS ${LAPACK_LIBRARY_DIR}
      NO_DEFAULT_PATH REQUIRED)

# Load lapack-config.cmake file
elseif (LAPACK_DIR OR LAPACK_ROOT)
   find_package(LAPACK CONFIG
      PATHS ${LAPACK_DIR} ${LAPACK_ROOT}
      NO_DEFAULT_PATH REQUIRED)

# Use cmake FindLAPACK.cmake module
else ()
   find_package(LAPACK MODULE REQUIRED)

endif ()


# Library is specified by the user
if (BLAS_LIBRARIES)
   if (NOT EXISTS ${BLAS_LIBRARIES})
      message(FATAL_ERROR "Blas library do not exist: ${BLAS_LIBRARIES}")
   endif ()

# Find directly the library
elseif (BLAS_LIBRARY_DIR)
   find_library(BLAS_LIBRARIES
      NAME blas openblas
      PATHS ${BLAS_LIBRARY_DIR}
      NO_DEFAULT_PATH REQUIRED)

# Load blas-config.cmake file
elseif (BLAS_DIR OR BLAS_ROOT)
   find_package(BLAS CONFIG
      PATHS ${BLAS_DIR} ${BLAS_ROOT}
      NO_DEFAULT_PATH REQUIRED)

# Use cmake FindBLAS.cmake module
else ()
   find_package(BLAS MODULE REQUIRED)

endif ()

message(VERBOSE "*** lapack libraries: ${LAPACK_LIBRARIES}")
message(VERBOSE "*** blas libraries: ${BLAS_LIBRARIES}")

set(LAPACK_BLAS_LIBRARIES "${LAPACK_LIBRARIES};${BLAS_LIBRARIES}")
