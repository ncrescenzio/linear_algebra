file(GLOB SUBROUTINE_FOLDER mi15*/)

if (EXISTS ${SUBROUTINE_FOLDER})
   message(STATUS "-- MI15 (Flexible GMRES) sources found")
   file(GLOB HSL_SOURCES ${SUBROUTINE_FOLDER}/src/*.f )
else()
   message(WARNING "MI15 (Flexible GMRES) sources NOT FOUND. Dummy sources are used.")
   file(GLOB HSL_SOURCES dummy_mi15/*.f)
endif()

add_library(hsl STATIC ${HSL_SOURCES})

include(GNUInstallDirs)

install(TARGETS hsl
   EXPORT linalg-targets
   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
   COMPONENT library
)
