cmake_minimum_required (VERSION 3.9)

add_compile_options(-std=c++11)
add_compile_options(-O3)

set(ROOT ${PROJECT_SOURCE_DIR})
set(SRC "${ROOT}/src")
include_directories("${SRC}")
include_directories("${SRC}/numericTools")

find_package(Eigen3 3.1.0 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_library(Quadrature ${SRC}/numericTools/quadratureRules.cpp)
add_library(Elliptic ${SRC}/numericTools/ellipticIntegral.cpp)
add_library(Legendre ${SRC}/numericTools/legendrePoly.cpp)



add_executable(testQuadrature ${ROOT}/dev/testQuadrature.cpp)
target_link_libraries(testQuadrature 
PUBLIC 
Quadrature
${EIGEN3_LIBS}
)

add_executable(testElliptic ${ROOT}/dev/testElliptic.cpp)
target_link_libraries(testElliptic 
PUBLIC 
Elliptic
${EIGEN3_LIBS}
)

add_executable(testLegendre ${ROOT}/dev/testLegendre.cpp)
target_link_libraries(testLegendre 
PUBLIC 
Legendre
Elliptic
${EIGEN3_LIBS}
)