cmake_minimum_required (VERSION 3.9)

option (TEST_NUMERICTOOLS       "" ON)
option (TEST_BEM      "" ON)
option (TEST_DEBUG      "" OFF)
option (TEST_CURRENT       "" ON)

if (TEST_DEBUG)
set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")
endif (TEST_DEBUG)


FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_CXX_FOUND)
message("OPENMP FOUND")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

endif()



add_compile_options(-std=c++11 -O3 -Wall -fopenmp)

set(ROOT ${PROJECT_SOURCE_DIR})
set(SRC "${ROOT}/src")
set(SRCBEM "${ROOT}/src/bem2D")
include_directories("${SRC}")
#include_directories("${SRC}")
#include_directories("${SRC}/numericTools")
#include_directories("${SRC}/spline")

find_package(Eigen3 3.1.0 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_library(NumericTools
${SRC}/numericTools/quadratureRules.cpp
${SRC}/numericTools/ellipticIntegral.cpp
${SRC}/numericTools/legendrePoly.cpp
${SRC}/numericTools/lagrangeBasis.cpp
${SRC}/numericTools/geometry2D.cpp
)

add_library(IO
${SRCBEM}/io/ioEigen.cpp
${SRCBEM}/io/utilities.cpp
)

add_library(BEM 
${SRCBEM}/boundaryElement.cpp 
${SRCBEM}/element.cpp 
${SRC}/spline/quinticSpline.cpp
)


if (TEST_NUMERICTOOLS)

add_executable(testQuadrature ${ROOT}/dev/testQuadrature.cpp)
target_link_libraries(testQuadrature 
PUBLIC 
NumericTools
${EIGEN3_LIBS}
)

add_executable(testElliptic ${ROOT}/dev/testElliptic.cpp)
target_link_libraries(testElliptic 
PUBLIC 
NumericTools
${EIGEN3_LIBS}
)

add_executable(testLegendre ${ROOT}/dev/testLegendre.cpp)
target_link_libraries(testLegendre 
PUBLIC 
NumericTools
${EIGEN3_LIBS}
)
endif (TEST_NUMERICTOOLS)

if (TEST_BEM)

add_executable(testBEMIntegral ${ROOT}/dev/testBEMIntegral.cpp)
target_link_libraries(testBEMIntegral 
PUBLIC 
BEM
NumericTools
IO
${EIGEN3_LIBS}
)

add_executable(testBEMAssemble ${ROOT}/dev/testBEMAssemble.cpp)
target_link_libraries(testBEMAssemble 
PUBLIC 
BEM
NumericTools
IO
${EIGEN3_LIBS}
OpenMP::OpenMP_CXX
)

add_executable(testBEMMixed ${ROOT}/dev/testBEMMixed.cpp)
target_link_libraries(testBEMMixed 
PUBLIC 
BEM
NumericTools
IO
${EIGEN3_LIBS}
OpenMP::OpenMP_CXX
)

endif (TEST_BEM)

add_executable(testBEMCone ${ROOT}/dev/testBEMCone.cpp)
target_link_libraries(testBEMCone 
PUBLIC 
BEM
NumericTools
IO
${EIGEN3_LIBS}
OpenMP::OpenMP_CXX
)
