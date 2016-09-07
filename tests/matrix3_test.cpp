#include "test.hpp"

namespace frea {
	namespace test {
		template <class T>
		using MatrixD_3 = RMatrix<T>;
		using TypesD_3 = ToTestTypes_t<types::SMatrixRange_t<types::Float_t, 3,4>>;
		TYPED_TEST_CASE(MatrixD_3, TypesD_3);

		TYPED_TEST(MatrixD_3, Rotation) {
			USING(mat_t);
			USING(value_t);
			using quat_t = QuatT<value_t, TestFixture::align>;

			const auto angle = this->makeRadian();
			constexpr auto Th = ThresholdF<value_t>(0.1);
			{
				// RotationX(ang) -> Rotation({1,0,0}, ang)
				const auto m0 = mat_t::RotationX(angle),
							m1 = mat_t::RotationAxis({1,0,0}, angle);
				EXPECT_LE(AbsMax(mat_t(m0-m1)), Th);
			}
			{
				// RotationY(ang) -> Rotation({0,1,0}, ang)
				const auto m0 = mat_t::RotationY(angle),
							m1 = mat_t::RotationAxis({0,1,0}, angle);
				EXPECT_LE(AbsMax(mat_t(m0-m1)), Th);
			}
			{
				// RotationZ(ang) -> Rotation({0,0,1}, ang)
				const auto m0 = mat_t::RotationZ(angle),
							m1 = mat_t::RotationAxis({0,0,1}, angle);
				EXPECT_LE(AbsMax(mat_t(m0-m1)), Th);
			}
			{
				// Rotation(axis, ang) -> Quaternion::Rotation(axis, ang)と比較
				const auto axis = this->makeDir();
				const auto m0 = mat_t::RotationAxis(axis, angle),
							m1 = quat_t::Rotation(axis, angle).asMat33();
				EXPECT_LE(AbsMax(mat_t(m0-m1)), ThresholdF<value_t>(0.5));
			}
		}
	}
}
