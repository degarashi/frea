#include "test.hpp"
#include "lubee/compiler_macro.hpp"
#include "../angle.hpp"

namespace frea {
	namespace test {
		template <class T>
		using VectorD_2 = RVector<T>;
		using TypesD_2 = ToTestTypes_t<types::VectorRangeD_t<types::Float_t, 2>>;
		TYPED_TEST_CASE(VectorD_2, TypesD_2);

		TYPED_TEST(VectorD_2, Wrap_Equality) {
			USING(vec_t);
			USING(value_t);
			const auto mtf = this->mt().template getUniformF<value_t>();
			const auto dir = random::GenVecUnitN<vec_t>(mtf, 2, 0.8);
			const auto w0 = dir[0].asInternal(),
						w1 = dir[1].asInternal();
			// cw
			EXPECT_EQ(dir[0].cw(dir[1]), w0.cw(w1));
			// ccw
			EXPECT_EQ(dir[0].ccw(dir[1]), w0.ccw(w1));
		}
		template <class V, class A>
		V Rotate(const V& v, const A& ang) {
			const double a = RadD(ang).get();
			const auto s = std::sin(a),
						c = std::cos(a);
			return V(v.x*c + v.y*s,
					v.x*-s + v.y*c);
		}
		TYPED_TEST(VectorD_2, Clockwise) {
			USING(value_t);
			using angle_t = Radian<value_t>;
			const auto v0 = this->makeRVec({-1e3, 1e3});
			const auto r_ang = [&mt=this->mt()](const lubee::Range<value_t>& r) {
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
