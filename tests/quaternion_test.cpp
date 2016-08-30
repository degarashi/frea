#include "test.hpp"
#include "../quaternion.hpp"
#include "../matrix.hpp"
#include "../random/angle.hpp"

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
					_rd(mt().template getUniformF<value_t>(DefaultRange))
				{}
				value_t makeRF() {
					return _rd();
				}
				rad_t makeRadian() {
					return random::GenHalfAngle<rad_t>(_rd);
				}
				vec_t makeVec3() {
					return random::GenVec<vec_t>(_rd);
				}
				vec_t makeDir() {
					return random::GenVecUnit<vec_t>(_rd);
				}
				quat_t makeRQuat() {
					return random::GenQuat<quat_t>(_rd);
				}
		};
		template <class T>
		constexpr Range<typename Quaternion<T>::value_t> Quaternion<T>::DefaultRange;

		TYPED_TEST_CASE(Quaternion, types::QTypes);

		TYPED_TEST(Quaternion, ConvertMatrix) {
			USING(mat3_t);
			USING(value_t);
			USING(quat_t);
			USING(array33_t);
			USING(vec_t);

			const auto ang = this->makeRadian();
			const auto axis = this->makeDir();
			auto q = quat_t::Rotation(axis, ang);
			const auto m = mat3_t::RotationAxis(axis, ang);

			// クォータニオンでベクトルを変換した結果が行列のそれと一致するか
			auto v = this->makeVec3();
			const auto v0 = vec_t(v * q),
						v1 = vec_t(v * m);
			constexpr auto Th = ThresholdF<value_t>(0.7);
            EXPECT_LT(AbsMax(vec_t(v0-v1)), Th);
			const array33_t ar0(m);
			// クォータニオンを行列に変換した結果が一致するか
			{
				const array33_t ar1(q.asMat33());
				EXPECT_LT(AbsMax(ar0 - ar1), Th);
			}
			// Matrix -> Quaternion -> Matrix の順で変換して前と後で一致するか
			q = quat_t::FromMat(m);
			{
				const array33_t ar1(q.asMat33());
				EXPECT_LT(AbsMax(ar0 - ar1), Th);
			}
		}
		TYPED_TEST(Quaternion, Multiply) {
			USING(quat_t);
			USING(vec_t);
			USING(mat3_t);
			USING(rad_t);
			USING(value_t);
			USING(array33_t);

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
			ASSERT_LT(AbsMax(array33_t(q[2].asMat33()) - m[2]), ThresholdF<value_t>(0.8));
		}
		TYPED_TEST(Quaternion, Rotation) {
			USING(quat_t);
			USING(vec_t);
			USING(rad_t);
			USING(value_t);

			// getRight(), getUp(), getDir()が{1,0,0},{0,1,0},{0,0,1}を変換した結果と比較
			const rad_t ang(this->makeRadian());
			const vec_t axis = this->makeDir();
			const auto q = quat_t::Rotation(axis, ang);
			const auto m = q.asMat33();

			constexpr auto Th = ThresholdF<value_t>(0.1);
			EXPECT_LT(AbsMax(vec_t(vec_t(1,0,0)*m - q.getRight())), Th);
			EXPECT_LT(AbsMax(vec_t(vec_t(0,1,0)*m - q.getUp())), Th);
			EXPECT_LT(AbsMax(vec_t(vec_t(0,0,1)*m - q.getDir())), Th);
		}
		//! クォータニオンの線形補間テスト
		TYPED_TEST(Quaternion, SLerp) {
			USING(vec_t);
			USING(quat_t);
			USING(rad_t);
			USING(value_t);
			USING(mat3_t);

			const int div = 8;
			const value_t tdiv = 1.0/div;
			const auto axis = this->makeDir();
			const rad_t ang(this->makeRadian());
			const auto q0 = this->makeRQuat();
			auto q1 = quat_t::Rotation(axis, ang);
			q1 = q0 * q1;
			const mat3_t m0 = q0.asMat33();
			for(int i=0 ; i<div ; i++) {
				const value_t t = tdiv * i;
				const mat3_t m1 = m0 * mat3_t::RotationAxis(axis, ang*t);
				const auto q2 = q0.slerp(q1, t);

				const auto v = this->makeDir();
				const vec_t v0 = v * q2;
				const vec_t v1 = v * m1;

				constexpr auto Th = ThresholdF<value_t>(0.9);
				EXPECT_LT(AbsMax(vec_t(v0 - v1)), Th);
			}
		}
	}
}
