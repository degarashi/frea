#include "test.hpp"

namespace frea {
	namespace test {
		template <class T>
		using VectorD_3 = RVector<T>;
		using TypesD_3 = ToTestTypes_t<types::VectorRangeD_t<types::FReg_t, 3>>;
		TYPED_TEST_CASE(VectorD_3, TypesD_3);

		// planeDivide
		// flip

		// TYPED_TEST(VectorD_3, Cross) {
		// }
		// TYPED_TEST(VectorD_3, VerticalVector) {
		// }
	}
}
