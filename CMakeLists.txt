cmake_minimum_required (VERSION 3.9)

option (TEST_NUMERICTOOLS       "" ON)
option (TEST_BEM      "" ON)
option (TEST_DEBUG      "" ON)

if (TEST_DEBUG)
set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")
endif (TEST_DEBUG)


add_compile_options(-std=c++11 -O3 -Wall)

set(ROOT ${PROJECT_SOURCE_DIR})
set(SRC "${ROOT}/src")
include_directories("${SRC}")
include_directories("${SRC}/numericTools")
include_directories("${SRC}/spline")

find_package(Eigen3 3.1.0 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})


add_library(NumericTools
${SRC}/numericTools/quadratureRules.cpp
${SRC}/numericTools/ellipticIntegral.cpp
${SRC}/numericTools/legendrePoly.cpp
${SRC}/numericTools/lagrangeBasis.cpp
${SRC}/numericTools/geometry2D.cpp
)

#add_library(Spline ${SRC}/spline/quinticSpline.cpp)
#add_library(Element ${SRC}/boundaryElement )
add_library(BEM ${SRC}/boundaryElement.cpp ${SRC}/element.cpp ${SRC}/spline/quinticSpline.cpp)


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

add_executable(testBEM ${ROOT}/dev/testBEM.cpp)
target_link_libraries(testBEM 
PUBLIC 
BEM
NumericTools
${EIGEN3_LIBS}
)
