cmake_minimum_required(VERSION 3.5)

project(hello)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
    set(optimised "-O2")
    set(full "-g -fbounds-check -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic      -fbacktrace -Wextra -Wconversion -Wconversion-extra")

endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

#set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${full}")
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${optimised}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

#
# Compile.
#
file(GLOB_RECURSE sources  src/*.f90 ../../multigrid/2D/multigrid.f90)
add_executable(CCAPS ${sources})

