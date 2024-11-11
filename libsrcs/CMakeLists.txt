
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

FILE(GLOB F90SOURCES ${CMAKE_CURRENT_LIST_DIR}/*.f90)
add_library(linalg ${F90SOURCES})

target_include_directories(linalg INTERFACE
   $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
   $<INSTALL_INTERFACE:mod>
)

target_link_libraries(linalg PRIVATE
   globals agmg hsl arpack ${LAPACK_BLAS_LIBRARIES})

include(GNUInstallDirs)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
   DESTINATION mod
   COMPONENT library
)
install(TARGETS linalg
   EXPORT linalg-targets
   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
   COMPONENT library
)
