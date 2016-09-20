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
				constexpr static lubee::Range<value_t> DefaultRange{-1e3, 1e3};
				constexpr static lubee::Range<value_t> DefaultAngleRange{-Radian<value_t>::OneRotationAng/2, Radian<value_t>::OneRotationAng/2};
				using RD = decltype(std::declval<Random>().mt().template getUniformF<value_t>(DefaultRange));
				RD	_rd;

			public:
				Quaternion():
					_rd(mt().template getUniformF<value_t>(DefaultRange))
				{}
				value_t makeRF(const lubee::Range<value_t>& r) {
					return mt().template getUniform<value_t>(r);
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
		constexpr lubee::Range<typename Quaternion<T>::value_t> Quaternion<T>::DefaultRange;

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
				constexpr auto Th = lubee::ThresholdF<value_t>(0.1);
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
			constexpr auto Th = lubee::ThresholdF<value_t>(0.3);
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
			constexpr auto Th = lubee::ThresholdF<value_t>(0.7);
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
		#define CHECK(m0, m1) \
			{ \
				auto q0 = this->makeRQuat(), \
					q1 = q0.m1(); \
					q0.m0(); \
					EXPECT_EQ(q0, q1); \
			}
		TYPED_TEST(Quaternion, FunctionEquality_Method) {
			CHECK(conjugate, conjugation)
			CHECK(invert, inversion)
			CHECK(normalize, normalization)
			{
				// rotate
				const auto axis = this->makeDir();
				const auto ang = this->makeRadian();
				auto q0 = this->makeRQuat(),
					 q1 = q0.rotation(axis, ang);
				q0.rotate(axis, ang);
				EXPECT_EQ(q0, q1);
			}
		}
		#undef CHECK
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
			q[2] = q[1] * q[0];
			q[2].normalize();
			m[2] = m[0] * m[1];
			ASSERT_LT(AbsMax(array33_t(q[2].asMat33()) - m[2]), lubee::ThresholdF<value_t>(0.8));
		}
		TYPED_TEST(Quaternion, Rotation) {
			USING(quat_t);
			USING(vec_t);
			USING(value_t);
			USING(rad_t);

			// getRight(), getUp(), getDir()が{1,0,0},{0,1,0},{0,0,1}を変換した結果と比較
			rad_t ang;
			do {
				ang = this->makeRadian();
			} while(std::abs(ang.get()) < 0.01);
			const auto axis = this->makeDir();
			const auto q = quat_t::Rotation(axis, ang);
			const auto m = q.asMat33();

			constexpr auto Th = lubee::ThresholdF<value_t>(0.1);
			EXPECT_LT(AbsMax(vec_t(vec_t(1,0,0)*m - q.getRight())), Th);
			EXPECT_LT(AbsMax(vec_t(vec_t(0,1,0)*m - q.getUp())), Th);
			EXPECT_LT(AbsMax(vec_t(vec_t(0,0,1)*m - q.getDir())), Th);

			// getXAxis() == 行列の1列目
			EXPECT_LE(AbsMax(vec_t(q.getXAxis() - m.template getColumn<0>())), Th);
			// Invは1行目
			EXPECT_LE(AbsMax(vec_t(q.getXAxisInv() - m.template getRow<0>())), Th);
			// getYAxis() == 行列の2列目
			EXPECT_LE(AbsMax(vec_t(q.getYAxis() - m.template getColumn<1>())), Th);
			// Invは2行目
			EXPECT_LE(AbsMax(vec_t(q.getYAxisInv() - m.template getRow<1>())), Th);
			// getZAxis() == 行列の3列目
			EXPECT_LE(AbsMax(vec_t(q.getZAxis() - m.template getColumn<2>())), Th);
			// Invは3行目
			EXPECT_LE(AbsMax(vec_t(q.getZAxisInv() - m.template getRow<2>())), Th);

			// RotationX(ang) == Rotation({1,0,0}, ang)
			EXPECT_EQ(quat_t::RotationX(ang), quat_t::Rotation({1,0,0}, ang));
			// RotationY(ang) == Rotation({0,1,0}, ang)
			EXPECT_EQ(quat_t::RotationY(ang), quat_t::Rotation({0,1,0}, ang));
			// RotationZ(ang) == Rotation({0,0,1}, ang)
			EXPECT_EQ(quat_t::RotationZ(ang), quat_t::Rotation({0,0,1}, ang));

			#define CHECK(func, x,y,z) { \
				auto q0 = q, \
					 q1 = q; \
				q0.func(ang); \
				q1.rotate({x,y,z}, ang); \
				EXPECT_EQ(q0, q1); \
			}
			// rotateX(ang) == rotate({1,0,0}, ang)
			CHECK(rotateX, 1,0,0)
			// rotateY(ang) == rotate({0,1,0}, ang)
			CHECK(rotateY, 0,1,0)
			// rotateZ(ang) == rotate({0,0,1}, ang)
			CHECK(rotateZ, 0,0,1)
			#undef CHECK

			{
				// angle & axis
				auto ang0 = q.angle().get(),
					 ang1 = ang.get();
				const auto axis0 = q.getAxis(),
					  axis1 = axis;
				constexpr auto Th = lubee::ThresholdF<value_t>(0.8);
				if(ang0 * ang1 < 0) {
					EXPECT_NEAR(axis0.dot(axis1), -1, Th);
					ang1 *= -1;
				}
				EXPECT_NEAR(ang0, ang1, Th);
			}
			{
				// rotation(from, to)
				const vec_t dir0 = this->makeDir();
				vec_t dir1;
				do {
					dir1 = this->makeDir();
				} while(std::abs(dir0.dot(dir1)) > 0.8);

				const auto q = quat_t::Rotation(dir0, dir1);
				const vec_t vdir = dir0 * q;
				constexpr auto Th = lubee::ThresholdF<value_t>(0.8);
				// Fromベクトルを変換したらToベクトルになる
				EXPECT_LE(AbsMax(vec_t(dir1-vdir)), Th);
				// 回転軸に平行なベクトルを回転させても値が変わらない
				const vec_t vvert0 = dir0.cross(dir1).normalization(),
							vvert1 = vvert0 * q;
				EXPECT_LE(AbsMax(vec_t(vvert1-vvert0)), Th);
			}
		}
		TYPED_TEST(Quaternion, YawPitchRoll) {
			USING(quat_t);
			USING(rad_t);
			const rad_t pitch(this->makeRF(rad_t::HalfRotationRange)),
						yaw(this->makeRF(rad_t::OneRotationRange)),
						roll(this->makeRF(rad_t::OneRotationRange));
			// Roll, Pitch, Yawの順番
			quat_t q0 = quat_t::Identity();
			q0.rotateZ(roll);
			q0.rotateX(-pitch);
			q0.rotateY(yaw);
			const quat_t q1 = quat_t::RotationYPR(yaw, pitch, roll);
			EXPECT_EQ(q0, q1);
		}
		TYPED_TEST(Quaternion, LookAt) {
			USING(vec_t);
			USING(value_t);
			USING(quat_t);
			// Z,Y軸方向のベクトルをLookAt(dir,up)で算出したクォータニオンで回転させ、それがdir,upと同一か確認
			const auto dir = this->makeDir();
			constexpr auto Th = lubee::ThresholdF<value_t>(0.5);
			vec_t up;
			do {
				up = this->makeDir();
			} while(std::abs(dir.dot(up)) > 1-Th);
			const auto q = quat_t::LookAt(dir, up);

			vec_t z(0,0,1),
				  y(0,1,0);
			z *= q;
			y *= q;
			EXPECT_NEAR(dir.dot(z), 1, Th);
			EXPECT_GE(up.dot(y), 0);
		}
		TYPED_TEST(Quaternion, SetLookAt) {
			USING(quat_t);
			USING(vec_t);
			USING(value_t);
			// 適当に重複しない軸フラグを決める
			Axis::e flag[3] = {Axis::X, Axis::Y, Axis::Z};
			std::shuffle(flag, flag+3, this->mt().refMt());
			// 重複しない位置座標をランダム生成
			vec_t pos = this->makeVec3(),
				  at, baseVec;
			constexpr auto Th = lubee::ThresholdF<value_t>(0.7);
			do {
				at = this->makeVec3();
			} while(pos.distance(at) < Th);
			// (at-pos)方向と重複しない方向ベクトルを生成
			const vec_t tdir = (at-pos).normalization();
			do {
				baseVec = this->makeDir();
			} while(std::abs(tdir.dot(baseVec)) > 0.9);
			const auto q = quat_t::SetLookAt(flag[0], flag[1], baseVec, at, pos);
			const vec_t axis[3] = {{1,0,0}, {0,1,0}, {0,0,1}};
			{
				// TargetAxisに指定された方向ベクトルは(at-pos)の方向と同じ
				const auto v = axis[flag[0]] * q;
				EXPECT_GE(v.dot(tdir), 1-Th);
			}
			{
				const auto v = axis[flag[1]] * q;
				// BaseAxisは多少補正される事があっても内積がマイナスにはならない
				EXPECT_GE(v.dot(baseVec), 0);
			}

			// 軸フラグが同一の時、InvalidAxis例外を送出
			for(int i=0 ; i<Axis::_Num ; i++)
				EXPECT_THROW(quat_t::SetLookAt(static_cast<Axis::e>(i), static_cast<Axis::e>(i), baseVec, at, pos), InvalidAxis);
			// atとposが同位置の時にNoValidAxis例外が送出される
			EXPECT_THROW(quat_t::SetLookAt(flag[0], flag[1], baseVec, at, at), NoValidAxis);
			// baseVecの方向が(at-pos)と同じか正対する時にNoValidAxisが送出される
			EXPECT_THROW(quat_t::SetLookAt(flag[0], flag[1], tdir, at, pos), NoValidAxis);
		}
		//! スカラとの四則演算, クォータニオンとの加減算
		TYPED_TEST(Quaternion, Operation) {
			USING(quat_t);
			const auto q = this->makeRQuat();
			const auto s = this->makeRF({-1e2, 1e2});

			#define CHECK(op) { \
				quat_t q0 = q op s, \
						q1 = q; \
				for(auto& qv : q1) \
					qv op##= s; \
				EXPECT_EQ(q0, q1); \
			}
			CHECK(+)
			CHECK(-)
			CHECK(*)
			CHECK(/)
			#undef CHECK

			const auto qa = this->makeRQuat();
			#define CHECK(op) { \
				quat_t q0 = q, \
						q1 = q op qa; \
				for(int i=0 ; i<4 ; i++) \
					q0[i] op##= qa[i]; \
				EXPECT_EQ(q0, q1); \
			}
			CHECK(+)
			CHECK(-)
			#undef CHECK
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
			q1 = q1 * q0;
			const mat3_t m0 = q0.asMat33();
			for(int i=0 ; i<div ; i++) {
				const value_t t = tdiv * i;
				const mat3_t m1 = m0 * mat3_t::RotationAxis(axis, ang*t);
				const auto q2 = q0.slerp(q1, t);

				const auto v = this->makeDir();
				const vec_t v0 = v * q2;
				const vec_t v1 = v * m1;

				constexpr auto Th = lubee::ThresholdF<value_t>(0.9);
				EXPECT_LT(AbsMax(vec_t(v0 - v1)), Th);
			}
			{
				// scale
				const auto t = this->makeRF({0.0, 1.0});
				const auto q_0 = quat_t::Identity().slerp(q0, t),
							q_1 = q0.scale(t);
				EXPECT_EQ(q_0, q_1);
			}
		}
		TYPED_TEST(Quaternion, Serialization) {
			CheckSerialization(this->makeRQuat());
		}
	}
}
