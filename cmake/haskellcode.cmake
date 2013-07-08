
file(GLOB HASKELL_SOURCES "${PROJECT_SOURCE_DIR}/*hs")

# add_custom_command(
#   OUTPUT haskell/functionals.exe
#   COMMAND cabal configure --builddir=${PROJECT_BINARY_DIR}/cabal && cabal build --builddir=${PROJECT_BINARY_DIR}/cabal
#   COMMENT ""
#   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src/haskell
#   DEPENDS ${PROJECT_SOURCE_DIR}/src/haskell/deft.cabal ${HASKELL_SOURCES}
# )


FILE(GLOB_RECURSE HASKELLSOURCE src/haskell/*hs)
set(HSD ${PROJECT_BINARY_DIR}/haskell)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/haskell/make.depend
  COMMAND rm -rf ${HSD}/make.depend && mkdir -p ${HSD} && touch ${HSD}/make.depend && HASDIR=${HSD}/ make depend
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src/haskell
  DEPENDS ${HASKELLSOURCE} ${PROJECT_SOURCE_DIR}/src/haskell/find-deps.pl ${PROJECT_SOURCE_DIR}/src/haskell/Makefile
)
add_custom_command(
  OUTPUT haskell/functionals.exe
  COMMAND echo && HASDIR=${HSD}/ make
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src/haskell
  DEPENDS ${HASKELLSOURCE} ${CMAKE_CURRENT_BINARY_DIR}/haskell/make.depend ${PROJECT_SOURCE_DIR}/src/haskell/Makefile
)
add_custom_target(
  functionals.exe
  DEPENDS haskell/functionals.exe
)

macro (haskell_functionals)
  foreach(filename ${ARGN})
    add_custom_command(
      OUTPUT ${filename}
      COMMAND mkdir -p src && haskell/functionals.exe ${filename}
      DEPENDS haskell/functionals.exe
    )
  endforeach(filename)
  add_library(defthaskell STATIC ${ARGN})
  target_link_libraries(defthaskell deftgeneric)
endmacro (haskell_functionals)

macro (new_haskell_functionals)
  foreach(filename ${ARGN})
    add_custom_command(
      OUTPUT ${filename}
      COMMAND mkdir -p src/new && haskell/newfunctionals.exe ${filename}
      DEPENDS haskell/functionals.exe
    )
  endforeach(filename)
  add_library(deftnewhaskell STATIC ${ARGN})
  target_link_libraries(deftnewhaskell deftgeneric)
endmacro (new_haskell_functionals)

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/doc)
macro (doc_pdf)
  set(TEMP_DOC_PDF_TEMP)
  foreach(filename ${ARGN})
    add_custom_command(
      OUTPUT ${filename}
      COMMAND haskell/latex-functionals.exe ${filename}
      DEPENDS haskell/functionals.exe
    )
    set(TEMP_DOC_PDF_TEMP ${TEMP_DOC_PDF_TEMP}; ${filename})
  endforeach(filename)
  add_custom_target(pdf
    DEPENDS ${TEMP_DOC_PDF_TEMP}
  )
endmacro (doc_pdf)
