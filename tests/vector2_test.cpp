#include "test.hpp"
#include "../compiler_macro.hpp"
#include "../angle.hpp"

namespace frea {
	namespace test {
		template <class T>
		using VectorD_2 = RVector<T>;
		using TypesD_2 = ToTestTypes_t<types::VectorRangeD_t<types::FReg_t, 2>>;
		TYPED_TEST_CASE(VectorD_2, TypesD_2);

		template <class V, class A>
		V Rotate(const V& v, const A& ang) {
			const double a = RadD(ang).get();
			const auto s = std::sin(a),
						c = std::cos(a);
			return V(v.x*c + v.y*s,
					v.x*-s + v.y*c);
		}
		TYPED_TEST(VectorD_2, Clockwise) {
			using value_t = typename TestFixture::value_t;
			using angle_t = Radian<value_t>;
			const auto v0 = this->makeRVec({-1e3, 1e3});
			const auto r_ang = [&mt=this->mt()](const Range<value_t>& r) {
				return mt.template getUniform<value_t>(r);
			};
			{
				// Degree(0, 180)ならClockwise
				const auto ang = r_ang({0, angle_t::OneRotationAng/2 * 0.99});
				const auto v1 = Rotate(v0, angle_t(ang));
				ASSERT_GE(v0.cw(v1), 0);
			}
			{
				// Degree(0, -180)ならCounterClockwise
				const auto ang = r_ang({-angle_t::OneRotationAng/2 * 0.99, 0});
				const auto v1 = Rotate(v0, angle_t(ang));
				ASSERT_GE(v0.ccw(v1), 0);
			}
			{
				// CcwとCwの結果は常に反対
				const auto ang = r_ang({-angle_t::OneRotationAng, angle_t::OneRotationAng});
				const auto v1 = Rotate(v0, angle_t(ang));
				ASSERT_LE(v0.cw(v1) * v0.ccw(v1), 0);
			}
		}
	}
}
