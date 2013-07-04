FILE(GLOB_RECURSE HASKELLSOURCE src/haskell/*hs)
add_custom_command(
  OUTPUT src/haskell/functionals.exe
  COMMAND make
  WORKING_DIRECTORY src/haskell
  DEPENDS ${HASKELLSOURCE}
)
add_custom_target(
  functionals.exe
  DEPENDS src/haskell/functionals.exe
)

macro (haskell_functionals)
  foreach(filename ${ARGN})
    add_custom_command(
      OUTPUT ${filename}
      COMMAND src/haskell/functionals.exe ${filename}
      DEPENDS src/haskell/functionals.exe
    )
  endforeach(filename)
  add_library(defthaskell STATIC ${ARGN})
  target_link_libraries(defthaskell deftgeneric)
endmacro (haskell_functionals)

macro (new_haskell_functionals)
  foreach(filename ${ARGN})
    add_custom_command(
      OUTPUT ${filename}
      COMMAND src/haskell/newfunctionals.exe ${filename}
      DEPENDS src/haskell/functionals.exe
    )
  endforeach(filename)
  add_library(deftnewhaskell STATIC ${ARGN})
  target_link_libraries(deftnewhaskell deftgeneric)
endmacro (new_haskell_functionals)
