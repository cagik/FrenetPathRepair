cmake_minimum_required(VERSION 3.0.2)
project(pc_test)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

add_definitions(-std=c++11 -w)

find_package(Eigen3 REQUIRED)

include_directories(
    include
    ${Eigen3_INCLUDE_DIRS}
    )

add_library(${PROJECT_NAME}
    src/KDTree.cpp
    src/mathfun.cpp
    src/pathrepairplanner.cpp
)

add_executable(${PROJECT_NAME}_demo src/main.cpp)
target_link_libraries(${PROJECT_NAME} /usr/local/lib/libosqp.so /usr/lib/x86_64-linux-gnu/libjsoncpp.so)
target_link_libraries(${PROJECT_NAME}_demo ${PROJECT_NAME})
