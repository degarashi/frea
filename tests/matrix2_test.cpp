#include "test.hpp"

namespace frea {
	namespace test {
		template <class T>
		using MatrixD_2 = RMatrix<T>;
		using TypesD_2 = ToTestTypes_t<types::SMatrixRange_t<types::Float_t, 2,3>>;
		TYPED_TEST_CASE(MatrixD_2, TypesD_2);

		TYPED_TEST(MatrixD_2, Rotation) {
			USING(mat_t);
			USING(vec_t);
			USING(value_t);

			constexpr auto range = Range<value_t>({-1e2, 1e2});
			const auto angle = this->makeRadian();
			const auto m = mat_t::Rotation(angle);
			const auto v0 = this->makeRVec(range);
			const auto fnRot = [](const auto& v, const auto& angv){
				return vec_t(v.x*std::cos(angv) - v.y*std::sin(angv),
							v.x*std::sin(angv) + v.y*std::cos(angv));
			};
			const vec_t v1a = v0 * m,
						v1b = fnRot(v0, angle.get());
			EXPECT_LE(AbsMax(vec_t(v1a-v1b)), ThresholdF<value_t>(0.5));
		}
	}
}
