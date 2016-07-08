#include "test.hpp"
#include "../compiler_macro.hpp"

namespace frea {
	namespace test {
		template <class T>
		using VectorD_4 = RVector<T>;
		using TypesD_4 = ToTestTypes_t<TypesD_t<4>>;
		TYPED_TEST_CASE(VectorD_4, TypesD_4);

		// asVec3Coord
	}
}
