#include <gtest/gtest.h>
#include "../matrix.hpp"
#include "../quaternion.hpp"
#include "../expquat.hpp"

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
