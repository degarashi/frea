#include "test.hpp"
#include "../quaternion.hpp"
#include "../matrix.hpp"

namespace frea {
	namespace test {
		using Types = ::testing::Types<std::tuple<float,BConst<false>>>;//, std::tuple<double, BConst<false>>>;
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
				std::function<value_t ()>	_rdf;

			public:
				Quaternion():
					_rdf(mt().template getUniformF<value_t>({-1e3, 1e3}))
				{}
				value_t makeRF() {
					return _rdf();
				}
				rad_t makeRadian() {
					return rad_t(makeRF());
				}
				vec_t makeVec3() {
					return random::GenVec<vec_t>(_rdf);
				}
				quat_t makeRQuat() {
					return random::GenQuat<quat_t>(_rdf);
				}
		};
		TYPED_TEST_CASE(Quaternion, Types);

		TYPED_TEST(Quaternion, ConvertMatrix) {
			using mat3_t = typename TestFixture::mat3_t;
			using value_t = typename TestFixture::value_t;
			using quat_t = typename TestFixture::quat_t;
			using array_t = typename TestFixture::array_t;
			using array33_t = typename TestFixture::array33_t;
			using vec_t = typename TestFixture::vec_t;

			const auto ang = this->makeRadian();
			const auto axis = this->makeVec3().normalization();
			auto q = quat_t::Rotation(axis, ang);
			const auto m = mat3_t::RotationAxis(axis, ang);

			// クォータニオンでベクトルを変換した結果が行列のそれと一致するか
			auto v = this->makeVec3();
			const array_t v0 = vec_t(v * q),
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
		// TYPED_TEST(Quaternion, Multiply) {
		// 	using This_t = std::decay_t<decltype(*this)>;
		// 	using QT = typename This_t::QuatType;
		// 	auto& rd = this->refRand();
		// 	auto rdf = rd.template getUniformF<float>();
		//
		// 	// クォータニオンを合成した結果を行列のケースと比較
		// 	for(int i=0 ; i<N_Iteration ; i++) {
		// 		RadF ang[2] = {RadF::Random(rdf), RadF::Random(rdf)};
		// 		Vec3 axis[2] = {Vec3::RandomDir(rdf), Vec3::RandomDir(rdf)};
		// 		QT		q[3];
		// 		AMat33	m[3];
		// 		for(int i=0 ; i<2 ; i++) {
		// 			q[i] = QT::Rotation(axis[i], ang[i]);
		// 			m[i] = AMat33::RotationAxis(axis[i], ang[i]);
		// 		}
		// 		q[2] = q[1] * q[0];
		// 		q[2].normalize();
		// 		m[2] = m[0] * m[1];
		// 		EXPECT_TRUE(Compare(q[2].asMat33(), m[2]));
		// 		// operator >> は * と逆の意味
		// 		q[2] = q[0] >> q[1];
		// 		q[2].normalize();
		// 		EXPECT_TRUE(Compare(q[2].asMat33(), m[2]));
		// 	}
		// }
		// TYPED_TEST(Quaternion, Rotation) {
		// 	using This_t = std::decay_t<decltype(*this)>;
		// 	constexpr bool Align = This_t::Align;
		// 	using QT = typename This_t::QuatType;
		// 	auto& rd = this->refRand();
		// 	auto rdf = rd.template getUniformF<float>();
		//
		// 	// getRight(), getUp(), getDir()が{1,0,0},{0,1,0},{0,0,1}を変換した結果と比較
		// 	for(int i=0 ; i<N_Iteration ; i++) {
		// 		RadF ang = RadF::Random(rdf);
		// 		Vec3 axis = Vec3::RandomDir(rdf);
		// 		QT q = QT::Rotation(axis, ang);
		// 		auto m = q.asMat33();
		//
		// 		EXPECT_TRUE(EqULPs(VecT<3,Align>(Vec3{1,0,0}*m), q.getRight(), ThresholdULPs));
		// 		EXPECT_TRUE(EqULPs(VecT<3,Align>(AVec3{0,1,0}*m), q.getUp(), ThresholdULPs));
		// 		EXPECT_TRUE(EqULPs(VecT<3,Align>(AVec3{0,0,1}*m), q.getDir(), ThresholdULPs));
		// 	}
		// }
		// TYPED_TEST(Quaternion, SLerp) {
		// 	using This_t = std::decay_t<decltype(*this)>;
		// 	using QT = typename This_t::QuatType;
		// 	auto& rd = this->refRand();
		//
		// 	// クォータニオンの線形補間
		// 	for(int i=0 ; i<N_Iteration ; i++) {
		// 		const int div = 32;
		// 		float tdiv = 1.f/div;
		// 		auto axis = this->genRandDir();
		// 		RadF ang = RadF::Random(rd.template getUniformF<float>());
		// 		auto q0 = this->genRandQ();
		// 		auto q1 = QT::Rotation(axis, ang);
		// 		q1 = q0 >> q1;
		// 		Mat33 m0 = q0.asMat33();
		// 		for(int i=0 ; i<div ; i++) {
		// 			float t = tdiv * i;
		// 			Mat33 m1 = m0 * Mat33::RotationAxis(axis, ang*t);
		// 			auto q2 = q0.slerp(q1, t);
		//
		// 			auto v = this->genRandDir();
		// 			auto v0 = v * q2;
		// 			auto v1 = v * m1;
		// 			EXPECT_TRUE(EqULPs(v0, v1, ThresholdULPs_Quat));
		// 		}
		// 	}
		// }
	}
}
