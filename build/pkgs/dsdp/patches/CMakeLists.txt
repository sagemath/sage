cmake_minimum_required(VERSION 3.0)
project (dsdp VERSION 5.8 LANGUAGES C)

link_directories(${LIBRARY_PREFIX}/lib)

find_package(LAPACK REQUIRED)

include_directories(include src/vecmat src/solver src/sdp)
FILE (GLOB_RECURSE   SRCS   src/*.c)
if(WIN32)
  list(APPEND SRCS dsdp.def)
endif()

## DSDP shared library
add_library(dsdp ${SRCS})
target_link_libraries(dsdp ${LAPACK_LIBRARIES})
target_include_directories(dsdp PUBLIC include src/vecmat src/solver src/sdp)
target_compile_options(dsdp PUBLIC -O2)

if(WIN32)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  target_compile_definitions(dsdp PUBLIC DSDP_MS_TIME)
else()
  target_link_libraries(dsdp m)
  target_compile_definitions(dsdp PUBLIC DSDP_TIME)
endif()

## DSDP5 executable
if(NOT WIN32)
  add_executable(dsdp5 examples/readsdpa.c)
  add_dependencies(dsdp5 dsdp)
  target_link_libraries(dsdp5 dsdp)
  target_include_directories(dsdp5 PUBLIC include)
endif()
if(APPLE)
  target_compile_definitions(dsdp5 PUBLIC _FORTIFY_SOURCE=0)
endif()

if(WIN32)
install(TARGETS dsdp
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
else()
install(TARGETS dsdp5 dsdp
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
endif()

FILE (GLOB DSDP_HEADERS include/*.h)
INSTALL(FILES ${DSDP_HEADERS} DESTINATION include)
