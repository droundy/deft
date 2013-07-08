enable_testing()

set(PAPER_DATS "")

file(GLOB PAPERS_SUBDIRECTOREIS papers/*/figs)
foreach(mypath ${PAPERS_SUBDIRECTOREIS})
  string(REPLACE ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} mygoodpath ${mypath})
  file(MAKE_DIRECTORY ${mygoodpath})
endforeach(mypath)

set(EXTRADATS)
# file(GLOB PAPERS_SUBDIRECTOREIS papers/*/figs/*.dat)
# foreach(mypath ${PAPERS_SUBDIRECTOREIS})
#   string(REPLACE ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} mygoodpath ${mypath})
#   add_custom_command(
#     OUTPUT ${mygoodpath}
#     COMMAND cp -v ${mypath} ${mygoodpath}
#     DEPENDS ${mypath}
#   )
#   set(EXTRADATS "${EXTRADATS}; ${mypath}")
# endforeach(mypath)

# message(STATUS "extradats is ${EXTRADATS}")

macro (add_papers_dat_using library)
  foreach(name ${ARGN})
    add_executable(papers/${name}.mkdat papers/${name}.cpp)
    target_link_libraries(papers/${name}.mkdat ${library})
    add_custom_command(
      OUTPUT papers/${name}.dat
      COMMAND papers/${name}.mkdat
      DEPENDS "papers/${name}.mkdat; ${EXTRADATS}"
    )
    set(PAPER_DATS ${PAPER_DATS}; papers/${name}.dat)
  endforeach(name)
endmacro (add_papers_dat_using)

macro (add_papers_mkdat_using library)
  foreach(name ${ARGN})
    add_executable(papers/${name}.mkdat papers/${name}.cpp)
    target_link_libraries(papers/${name}.mkdat ${library})
    #add_custom_command(
    #  OUTPUT papers/${name}.dat
    #  COMMAND papers/${name}.mkdat
    #  DEPENDS papers/${name}.mkdat
    #)
  endforeach(name)
endmacro (add_papers_mkdat_using)
