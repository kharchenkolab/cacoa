cmake_minimum_required(VERSION 2.8)

project(cacoa)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -std=c++11 -Wno-attributes -Wno-deprecated -Wno-deprecated-declarations")
set(FIND_RCPP_PROGRESS ON)

find_package (R REQUIRED)

set(INCLUDE_DIRS ${R_INCLUDE_DIRS})
include_directories(${INCLUDE_DIRS})

set(BASE_LIBRARIES ${R_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

FILE(GLOB CacoaSources ./*.cpp)
add_library(Cacoa ${CacoaSources})
target_link_libraries(Cacoa ${BASE_LIBRARIES})
