cmake_minimum_required(VERSION 2.8)
project(frea)
include(aux.cmake)

set(CXX_VERSION 17)
include(lubee/cxx_options.cmake)
string(CONCAT CMAKE_CXX_FLAGS " -Wno-gnu-anonymous-struct")
string(CONCAT CMAKE_CXX_FLAGS " -Wno-nested-anon-types")

add_definitions(-DSSE=${SSE})
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
)

if(NOT without-test)
	enable_testing()
	find_package(Threads REQUIRED)
	find_package(GTest REQUIRED)
	include_directories(${GTEST_INCLUDE_DIRS})
	# testsディレクトリ以下のソースが対象
	aux_source_directory(src/tests TEST_SRC)
	# matrix, vector関連の物は個別にexeを作成
	MakeSeparationTest("matrix;vector" TEST_SRC)
	# 他はotherとしてまとめる
	AddTest(other "${TEST_SRC}")
endif()
