cmake_minimum_required(VERSION 2.8)
project(frea)

find_package(Threads REQUIRED)
find_package(GTest REQUIRED)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -ggdb3 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_DEBUG} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++14 -Wall -Wextra -ftemplate-depth=1024")
add_definitions(-DSSE=${SSE})
include_directories(${GTEST_INCLUDE_DIRS})
aux_source_directory(tests TEST_SRC)
add_executable(frea_test ${TEST_SRC})
target_link_libraries(frea_test ${CMAKE_THREAD_LIBS_INIT}
								${GTEST_LIBRARIES})
