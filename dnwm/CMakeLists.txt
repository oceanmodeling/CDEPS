project(dnwm Fortran)
set(SRCFILES  nwm_comp_nuopc.F90)

foreach(FILE ${SRCFILES})
  if(EXISTS "${CASEROOT}/SourceMods/src.dnwm/${FILE}")
    list(REMOVE_ITEM SRCFILES ${FILE})
    list(APPEND SRCFILES "${CASEROOT}/SourceMods/src.dnwm/${FILE}")
    message("Using ${FILE} from ${CASEROOT}/SourceMods/src.dnwm")
  endif()
endforeach()

message("DNWM srcfiles are ${SRCFILES}")

add_library(dnwm ${SRCFILES})

add_dependencies(dnwm dshr streams)
target_include_directories (dnwm PRIVATE ${ESMF_F90COMPILEPATHS})
target_include_directories (dnwm PRIVATE "${CMAKE_SOURCE_DIR}")
target_include_directories (dnwm PRIVATE "${PIO_Fortran_INCLUDE_DIR}")
if(NOT DISABLE_FoX)
target_include_directories (dnwm PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../fox/include)
endif()

if(BLD_STANDALONE)
  # ESMX requires mod files
  foreach (SRC ${SRCFILES})
    string(REGEX REPLACE "[.]F90$" ".mod" MOD ${SRC})
    if (NOT DEFINED CIMEROOT AND MOD STREQUAL nwm_comp_nuopc.mod)
      set(MOD cdeps_dnwm_comp.mod)
    endif()
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${MOD}" DESTINATION include)
  endforeach ()
endif()
