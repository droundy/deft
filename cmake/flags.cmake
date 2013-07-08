# to adjust flags nicely:

# cmake -DCMAKE_BUILD_TYPE=Release .
# cmake -DCMAKE_BUILD_TYPE=Debug .

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ansi -W -Wall -Werror -O2 -Wno-unused-variable -Wno-unused-parameter -Wno-return-type -Wno-unused-local-typedefs")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ansi")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY libs)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY libs)
