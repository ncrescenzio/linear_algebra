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
