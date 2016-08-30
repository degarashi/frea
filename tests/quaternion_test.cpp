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
				using array44_t = ArrayM<value_t, 4, 4>;
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

		TYPED_TEST(Quaternion, Identity) {
			USING(quat_t);
			USING(value_t);
			{
				// 単位クォータニオンをどちらから掛けても値が変わらない
				const auto q0 = this->makeRQuat(),
							qi = quat_t::Identity();
				const auto q1 = qi * q0,
							q2 = q0 * qi;
				constexpr auto Th = ThresholdF<value_t>(0.1);
				ASSERT_LE(q0.distance(q1), Th);
				ASSERT_LE(q0.distance(q2), Th);
			}
		}
		TYPED_TEST(Quaternion, Compare) {
			USING(value_t);

			const auto q0 = this->makeRQuat();
			auto q1 = q0;
			ASSERT_EQ(q0, q0);
			ASSERT_EQ(q0, q1);
			auto& mt = this->mt();
			q1[mt.template getUniform<int>({0,3})] += mt.template getUniform<value_t>();
			ASSERT_NE(q0, q1);
		}
		TYPED_TEST(Quaternion, Invert) {
			USING(value_t);
			USING(quat_t);
			// 逆クォータニオンをどちらから掛けても単位クォータニオンになる
			const auto q0 = this->makeRQuat(),
						q1 = q0.inversion(),
						qi = quat_t::Identity();
			constexpr auto Th = ThresholdF<value_t>(0.3);
			ASSERT_LE((q0*q1).distance(qi), Th);
			ASSERT_LE((q1*q0).distance(qi), Th);
		}
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
			{
				USING(array44_t);
				USING(mat4_t);
				// 4x4行列のテスト
				const array44_t ar1(q.asMat44()),
								ar0(mat4_t(m.template convertI<4,4,3>(1)));
				EXPECT_LT(AbsMax(ar0 - ar1), Th);
			}
			// Matrix -> Quaternion -> Matrix の順で変換して前と後で一致するか
			q = quat_t::FromMat(m);
			{
				const array33_t ar1(q.asMat33());
				EXPECT_LT(AbsMax(ar0 - ar1), Th);
			}
		}
		TYPED_TEST(Quaternion, FunctionEquality_Static) {
			USING(quat_t);
			{
				// Identity
				quat_t q0, q1=quat_t::Identity();
				q0.identity();
				EXPECT_EQ(q0, q1);
			}
			const auto angle = this->makeRadian();
			{
				// RotationXYZ -> rotationXYZ
				const quat_t q0[3] = {
					quat_t::RotationX(angle),
					quat_t::RotationY(angle),
					quat_t::RotationZ(angle)
				};
				quat_t q1[3];
				for(auto& q : q1)
					q = quat_t::Identity();
				q1[0] = q1[0].rotationX(angle);
				q1[1] = q1[1].rotationY(angle);
				q1[2] = q1[2].rotationZ(angle);

				for(int i=0 ; i<3 ; i++)
					EXPECT_EQ(q0[i], q1[i]);
			}
			const auto axis = this->makeDir();
			{
				// Rotation -> rotation
				const auto q0 = quat_t::Rotation(axis, angle);
				auto q1 = quat_t::Identity();
				q1.rotate(axis, angle);
				EXPECT_EQ(q0, q1);
			}
		}
		TYPED_TEST(Quaternion, FunctionEquality_Method) {
			USING(quat_t);
			{
				// conjugate
				auto q0 = this->makeRQuat(),
					 q1 = q0.conjugation();
				q0.conjugate();
				EXPECT_EQ(q0, q1);
			}
			{
				// invert
				auto q0 = this->makeRQuat(),
					 q1 = q0.inversion();
				q0.invert();
				EXPECT_EQ(q0, q1);
			}
			{
				// normalize
				auto q0 = this->makeRQuat(),
					 q1 = q0.normalization();
				q0.normalize();
				EXPECT_EQ(q0, q1);
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
