cmake_minimum_required(VERSION 2.8)
project(frea)

set(COMMON_COMPILE lubee/common_compile)
include(${COMMON_COMPILE}/cxx_setup.cmake)
include(${COMMON_COMPILE}/compile_commands.cmake)
include(${COMMON_COMPILE}/aux.cmake)
SetupCXX(17)
string(CONCAT CMAKE_CXX_FLAGS " -Wno-gnu-anonymous-struct")
string(CONCAT CMAKE_CXX_FLAGS " -Wno-nested-anon-types")

add_definitions(-DSSE=${SSE})
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
)

if(NOT without-test)
	enable_testing()
	# testsディレクトリ以下のソースが対象
	aux_source_directory(src/tests TEST_SRC)
	LoadGTestLibs(LIBS)
	# matrix, vector関連の物は個別にexeを作成
	MakeSeparationTest(${PROJECT_NAME} "matrix;vector" TEST_SRC "${LIBS}")
	# 他はotherとしてまとめる
	AddTest(${PROJECT_NAME} other "${TEST_SRC}" "${LIBS}")
	DefineCompDB(TEST_SRC)
endif()
