file(GLOB local arpack-ng/SRC)

if (EXISTS ${local})
   message(STATUS "-- Arpack FOUND.")
   file(GLOB ARPACK_SOURCES arpack-ng/SRC/*.f arpack-ng/UTIL/*.f)
   add_library(arpack STATIC ${ARPACK_SOURCES})
else()
   message(WARNING "Arpack sources NOT FOUND. Dummy sources are used.")
   add_library(arpack STATIC dummy_arpack.f)
endif()

install(TARGETS arpack
   EXPORT linalg-targets
   ARCHIVE DESTINATION lib
   COMPONENT library
)
