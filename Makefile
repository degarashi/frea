PWD			:= $(shell pwd)
WORK_DIR	:= /var/tmp/frea

CMake = mkdir -p $(1); cd $(1); cmake $(PWD) -G 'CodeBlocks - Unix Makefiles' -DSSE=2 -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=clang++;
Make = cd $(1); make -j$(2);
Clean = cd $(1); make clean;
all:
	$(call CMake, $(WORK_DIR))
	$(call Make, $(WORK_DIR), 5)
clean:
	$(call Clean, $(WORK_DIR))
