file(GLOB local AGMG*/)

if (EXISTS ${local})
   message(STATUS "-- AGMG sources FOUND")
   file(GLOB agmg_SOURCES ${local}/SRC/dagmg.f90 ${local}/SRC/dagmg_mumps.f90)
else()
   message(WARNING "AGMG sources NOT FOUND. Dummy sources are used instead. Library will not include AGMG tools ")
   file(GLOB agmg_SOURCES dummy_src/*.f90)
endif()

add_library(agmg STATIC ${agmg_SOURCES})
