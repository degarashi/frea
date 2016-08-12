#include "test.hpp"
#include "../quaternion.hpp"
#include "../matrix.hpp"

namespace frea {
	namespace test {
		template <class T>
		class Quaternion : public Random {
			public:
				using value_t = typename std::tuple_element_t<0,T>;
				constexpr static bool align = std::tuple_element_t<1,T>::value;
				using quat_t = QuatT<value_t, align>;
				using rad_t = typename quat_t::rad_t;
				using vec_t = typename quat_t::vec_t;
				using vec4_t = typename quat_t::vec4_t;
				using mat3_t = typename quat_t::mat3_t;
				using mat4_t = typename quat_t::mat4_t;
				using exp_t = typename quat_t::exp_t;
				using array_t = Array<value_t, 3>;
				using array33_t = ArrayM<value_t, 3, 3>;
			private:
				constexpr static Range<value_t> DefaultRange{-1e3, 1e3};
				constexpr static Range<value_t> DefaultAngleRange{-Radian<value_t>::OneRotationAng/2, Radian<value_t>::OneRotationAng/2};
				using RD = decltype(std::declval<Random>().mt().template getUniformF<value_t>(DefaultRange));
				RD	_rd;

			public:
				Quaternion():
					_rd(mt().template getUniformF<value_t>())
				{}
				value_t makeRF() {
					constexpr auto R = DefaultRange;
					return mt().template getUniform<value_t>(R);
				}
				rad_t makeRadian() {
					constexpr auto R = DefaultAngleRange;
					return rad_t(mt().template getUniform<value_t>(R));
				}
				vec_t makeVec3() {
					return random::GenVec<vec_t>(_rd);
				}
				vec_t makeDir() {
					return random::GenDir<vec_t>(_rd);
				}
				quat_t makeRQuat() {
					return random::GenQuat<quat_t>(_rd);
				}
		};
		TYPED_TEST_CASE(Quaternion, types::QTypes);

		TYPED_TEST(Quaternion, ConvertMatrix) {
			using mat3_t = typename TestFixture::mat3_t;
			using value_t = typename TestFixture::value_t;
			using quat_t = typename TestFixture::quat_t;
			using array33_t = typename TestFixture::array33_t;
			using vec_t = typename TestFixture::vec_t;

			const auto ang = this->makeRadian();
			const auto axis = this->makeVec3().normalization();
			auto q = quat_t::Rotation(axis, ang);
			const auto m = mat3_t::RotationAxis(axis, ang);

			// クォータニオンでベクトルを変換した結果が行列のそれと一致するか
			auto v = this->makeVec3();
			const auto v0 = vec_t(v * q),
						v1 = vec_t(v * m);
			EXPECT_TRUE(ulps::Equal(v0, v1, ThresholdULPs<value_t>));
			const array33_t ar0(m);
			// クォータニオンを行列に変換した結果が一致するか
			{
				const array33_t ar1(q.asMat33());
				EXPECT_LT(AbsMax(ar0 - ar1), 1e-2);
			}
			// Matrix -> Quaternion -> Matrix の順で変換して前と後で一致するか
			q = quat_t::FromMat(m);
			{
				const array33_t ar1(q.asMat33());
				EXPECT_LT(AbsMax(ar0 - ar1), 1e-2);
			}
		}
		TYPED_TEST(Quaternion, Multiply) {
			using quat_t = typename TestFixture::quat_t;
			using vec_t = typename TestFixture::vec_t;
			using mat3_t = typename TestFixture::mat3_t;
			using rad_t = typename TestFixture::rad_t;
			using array33_t = typename TestFixture::array33_t;

			// クォータニオンを合成した結果を行列のケースと比較
			const rad_t ang[2] = {this->makeRadian(), this->makeRadian()};
			const vec_t axis[2] = {this->makeDir(), this->makeDir()};
			quat_t	q[3];
			mat3_t	m[3];
			for(int i=0 ; i<2 ; i++) {
				q[i] = quat_t::Rotation(axis[i], ang[i]);
				m[i] = mat3_t::RotationAxis(axis[i], ang[i]);
			}
			q[2] = q[0] * q[1];
			q[2].normalize();
			m[2] = m[0] * m[1];
			ASSERT_LT(AbsMax(array33_t(q[2].asMat33()) - m[2]), 1e-2);
		}
		TYPED_TEST(Quaternion, Rotation) {
			using quat_t = typename TestFixture::quat_t;
			using vec_t = typename TestFixture::vec_t;
			using rad_t = typename TestFixture::rad_t;
			using value_t = typename TestFixture::value_t;

			// getRight(), getUp(), getDir()が{1,0,0},{0,1,0},{0,0,1}を変換した結果と比較
			const rad_t ang(this->makeRadian());
			const vec_t axis = this->makeDir();
			const auto q = quat_t::Rotation(axis, ang);
			const auto m = q.asMat33();

			EXPECT_TRUE(ulps::Equal(vec_t(vec_t(1,0,0)*m), q.getRight(), ThresholdULPs<value_t>));
			EXPECT_TRUE(ulps::Equal(vec_t(vec_t(0,1,0)*m), q.getUp(), ThresholdULPs<value_t>));
			EXPECT_TRUE(ulps::Equal(vec_t(vec_t(0,0,1)*m), q.getDir(), ThresholdULPs<value_t>));
		}
		//! クォータニオンの線形補間テスト
		TYPED_TEST(Quaternion, SLerp) {
			using vec_t = typename TestFixture::vec_t;
			using quat_t = typename TestFixture::quat_t;
			using rad_t = typename TestFixture::rad_t;
			using value_t = typename TestFixture::value_t;
			using mat3_t = typename TestFixture::mat3_t;
			constexpr auto ThresholdULPs_Quat = ulps::Diff_C<value_t>(0.0, 5e-3);

			const int div = 8;
			const value_t tdiv = 1.0/div;
			const auto axis = this->makeDir();
			const rad_t ang(this->makeRadian());
			const auto q0 = this->makeRQuat();
			auto q1 = quat_t::Rotation(axis, ang);
			q1 = q0 * q1;
			const mat3_t m0 = q0.asMat33();
			for(int i=0 ; i<div ; i++) {
				value_t t = tdiv * i;
				mat3_t m1 = m0 * mat3_t::RotationAxis(axis, ang*t);
				auto q2 = q0.slerp(q1, t);

				auto v = this->makeDir();
				const vec_t v0 = v * q2;
				const vec_t v1 = v * m1;
				EXPECT_TRUE(ulps::Equal(v0, v1, ThresholdULPs_Quat));
			}
		}
	}
}
