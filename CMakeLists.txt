# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project( reconstruct )
set(CMAKE_CXX_LINK_FLAGS "-lGL")
set(CMAKE_CXX_FLAGS, "-g -pthread -Ofast")
set(CMAKE_BUILD_TYPE "Release")

# SET YOUR VTK and CGAL directory here!!
SET(VTK_DIR /home/gnomble/VTK-9.1.0/build/lib/cmake/vtk-9.1)
SET(CGAL_DIR /home/gnomble/CGAL-5.3.1/lib/cmake/CGAL)

# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #

# CGAL and its components
find_package( CGAL REQUIRED  )
include( ${CGAL_USE_FILE})

# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #

# VTK 
find_package( VTK REQUIRED  )
#include( ${VTK_USE_FILE})

# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #

# Boost and its components
find_package( Boost REQUIRED )

# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #

include(CGAL_CreateSingleSourceCGALProgram)
include_directories(BEFORE "../../include")
add_executable(reconstruct reconstruct.cxx util.cxx)

# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #

#add_executable(reconstruct reconstruct.cxx)

# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #

target_link_libraries(reconstruct  ${VTK_LIBRARIES} )

# Creating entries for all C++ files with "main" routine
# ##########################################################


vtk_module_autoinit(
    TARGETS reconstruct 
    MODULES ${VTK_LIBRARIES}
    )


