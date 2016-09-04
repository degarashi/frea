#include "test.hpp"

namespace frea {
	namespace test {
		template <class T>
		using FMatrix = RMatrix<T>;
		using Types_t = types::MatrixRange_t<types::FReg_t, 3,5, 3,5>;
		using Types = ToTestTypes_t<Types_t>;
		TYPED_TEST_CASE(FMatrix, Types);

		TYPED_TEST(FMatrix, LinearNormalize) {
			USING(value_t);
			USING(mat_t);
			USING(vec_t);

			constexpr auto range = Range<value_t>{-1e2, 1e2};
			const auto m = this->makeRMat(range),
						mn0 = m.linearNormalization();
			value_t r[mat_t::dim_m];
			for(int i=0 ; i<mat_t::dim_m ; i++)
				r[i] = m.m[i].absolute().getMaxValue();
			// 割った数をかければ元と大体同じ(行ごとに比較する)
			for(int i=0 ; i<mat_t::dim_m ; i++) {
				const vec_t& v0 = m.m[i];
				const vec_t v1 = mn0.m[i] * r[i];
				EXPECT_LE(AbsMax(vec_t(v1-v0)), ThresholdF<value_t>(0.3));
			}
		}
	}
}
