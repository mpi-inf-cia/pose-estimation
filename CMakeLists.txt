cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
project(tracking)
include_directories(graph/include graphics/include include)
file(GLOB headers
    include/andres/*.hxx
    include/andres/graph/*.hxx
    include/andres/graph/hdf5/*.hxx
    include/andres/graphics/*.hxx
    include/andres/ilp/*.hxx
    include/pose/*.hxx
)
enable_testing()

##############################################################################
# GNU g++ specific settings
##############################################################################
if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=c++11)
endif()

##############################################################################
# OpenMP
##############################################################################
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP found")
else()
    message("OpenMP not found")
endif()

##############################################################################
# HDF5
##############################################################################
find_package(HDF5)
if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
endif()

##############################################################################
# Gurobi
##############################################################################
find_package(GUROBI)
if(GUROBI_FOUND)
    include_directories(${GUROBI_INCLUDE_DIR})
endif()

##############################################################################
# GLUT
##############################################################################
find_package(GLUT)

##############################################################################
# targets: graphics
##############################################################################
if(HDF5_FOUND AND GLUT_FOUND)
    add_executable(view3d src/andres/graphics/viewer-opengl.cxx ${headers})
    target_link_libraries(view3d ${HDF5_LIBRARIES} ${GLUT_LIBRARIES})
endif()

##############################################################################
# targets: pose
##############################################################################
# test
add_executable(test-problem src/pose/test/problem.cxx ${headers})
add_test(test-problem test-problem)

add_executable(test-solver src/pose/test/solver.cxx ${headers})
add_test(test-solver test-solver)

# research
if(HDF5_FOUND)
    add_executable(plot-problem src/pose/research/plot-problem.cxx ${headers})
    target_link_libraries(plot-problem ${HDF5_LIBRARIES})

    if(GUROBI_FOUND)
        add_executable(solver src/pose/research/solver.cxx ${headers})
        target_link_libraries(solver ${HDF5_LIBRARIES} ${GUROBI_LIBRARIES})

        add_executable(solver-callback src/pose/research/solver-callback.cxx ${headers})
        target_link_libraries(solver-callback ${HDF5_LIBRARIES} ${GUROBI_LIBRARIES})
    endif()
endif()
