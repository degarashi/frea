#include "test.hpp"
#include "../expquat.hpp"

namespace frea {
	namespace test {
		template <class T>
		class ExpQuaternion : public Random {
			public:
				using value_t = typename std::tuple_element_t<0,T>;
				constexpr static bool align = std::tuple_element_t<1,T>::value;
				using eq_t = ExpQuatT<value_t, align>;
				using quat_t = QuatT<value_t, align>;
				using vec_t = typename quat_t::vec_t;
				using array_t = Array<value_t, 3>;

			private:
				constexpr static Range<value_t> DefaultRange{-1e3, 1e3};
				using RD = decltype(std::declval<Random>().mt().template getUniformF<value_t>(DefaultRange));
				RD	_rd;

			public:
				ExpQuaternion():
					_rd(mt().template getUniformF<value_t>())
				{}
				vec_t makeDir() {
					return random::GenDir<vec_t>(_rd);
				}
				quat_t makeRQuat() {
					return random::GenQuat<quat_t>(_rd);
				}
		};
		TYPED_TEST_CASE(ExpQuaternion, types::QTypes);

		TYPED_TEST(ExpQuaternion, QuaternionConvert) {
			using eq_t = typename TestFixture::eq_t;
			using vec_t = typename TestFixture::vec_t;
			using value_t = typename TestFixture::value_t;

			// Quat -> ExpQuat -> Quat の結果比較
			const auto q0 = this->makeRQuat();
			eq_t eq(q0);
			const auto q1 = eq.asQuat();
			eq = eq_t(q1*-1);
			const auto q2 = eq.asQuat();
			// クォータニオンの符号反転した物は同一視 = 方向ベクトルを変換して同じだったらOK
			const auto v0 = this->makeDir();
			const vec_t v1 = v0 * q0,
						v2 = v0 * q1,
						v3 = v0 * q2;
			constexpr value_t Threshold = 8e-2;
			EXPECT_LT(AbsMax(vec_t(v1-v2)), Threshold);
			EXPECT_LT(AbsMax(vec_t(v1-v3)), Threshold);
		}
		TYPED_TEST(ExpQuaternion, Lerp) {
			using eq_t = typename TestFixture::eq_t;
			using value_t = typename TestFixture::value_t;
			// ExpQuatを合成した結果をQuat合成のケースと比較
			auto q0 = this->makeRQuat(),
				q1 = this->makeRQuat();
			// 最短距離で補間するような細工
			if(q0.dot(q1) < 0)
				q1 *= -1;
			eq_t eq0(q0),
				eq1(q1),
				eq10(q0 * -1),
				eq11(q1 * -1);
			const auto len0 = eq0.asVec3().length(),
						len1 = eq1.asVec3().length(),
						len10 = eq10.asVec3().length(),
						len11 = eq11.asVec3().length();
			if(len0+len1 > len10+len11) {
				eq0 = q0 * -1;
				eq1 = q1 * -1;
				eq10 = q0;
				eq11 = q1;
			}

			constexpr value_t Threshold = 0.4;
			constexpr int NDiv = 8;
			value_t err_sum = 0;
			const auto v0 = this->makeDir();
			// 0から1まで遷移
			for(int i=0 ; i<=NDiv ; i++) {
				const value_t t = i / value_t(NDiv);
				// Quatで変換
				const auto v1 = v0 * q0.slerp(q1, t);
				// ExpQuatで変換
				const auto eq2 = eq0 * (1.0 - t) + eq1 * t;
				const auto v2 = v0 * eq2.asQuat();

				const auto diff = v2-v1;
				err_sum += diff.dot(diff);
				EXPECT_LT(err_sum/(i+1), Threshold);
			}
			err_sum /= NDiv+1;
			EXPECT_LT(err_sum, Threshold);
		}
	}
}
