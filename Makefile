LIB_NAME			:= frea
MAKE_GDBINIT		:= YES
SSE					?= 2

OPT_SSE					= -DSSE=$(SSE)

ADDITIONAL_CMAKE_OPTION	:= $(OPT_SSE)
include lubee/common_compile/common.make
