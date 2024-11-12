# Search "globals" library
if(NOT TARGET "globals")

   # Try to find "globals" library using
   # -Dglobals_DIR, -Dglobals_ROOT, -DGLOBALS_DIR and -DGLOBALS_ROOT
   find_package(globals PATHS ${GLOBALS_ROOT} ${GLOBALS_DIR} QUIET)

   if(NOT globals_FOUND)

      # Assume that source code for globals library
      # is in the same directory of linalg library
      set(GLOBALS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../globals)

      if(EXISTS ${GLOBALS_SRC})

         # Assume the code is already compiled
         if(EXISTS ${GLOBALS_SRC}/build)
            find_package(globals REQUIRED
               PATHS ${GLOBALS_SRC}/build/
               NO_DEFAULT_PATH)

         # The code is not compiled
         else(EXISTS ${GLOBALS_SRC}/build)
            add_subdirectory(${GLOBALS_SRC} ${GLOBALS_SRC}/build EXCLUDE_FROM_ALL)

         endif(EXISTS ${GLOBALS_SRC}/build)

      else(EXISTS ${GLOBALS_SRC})
         message(FATAL_ERROR "Library 'globals' NOT FOUND!")

      endif(EXISTS ${GLOBALS_SRC})

   endif(NOT globals_FOUND)

endif(NOT TARGET "globals")
