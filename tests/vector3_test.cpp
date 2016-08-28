#include "test.hpp"

namespace frea {
	namespace test {
		template <class T>
		using VectorD_3 = RVector<T>;
		using TypesD_3 = ToTestTypes_t<types::VectorRangeD_t<types::FReg_t, 3>>;
		TYPED_TEST_CASE(VectorD_3, TypesD_3);

		namespace {
			template <class value_t, class MT>
			value_t Scale(MT& mt) {
				const value_t s = mt.template getUniform<int>({0,1}) == 0 ? -1 : 1;
				return s * mt.template getUniform<value_t>({1e-1, 1e1});
			}
			template <class vec_t, class MT>
			auto GenVecPair(MT& mt) {
				using value_t = typename vec_t::value_t;
				// 同じ方向か、正反対のベクトルの組は避ける
				const auto mtf = mt.template getUniformF<value_t>();
				vec_t v0 = random::GenVecUnit<vec_t>(mtf),
					 v1;
				for(;;) {
					v1 = random::GenVecUnit<vec_t>(mtf);
					const auto d = v0.dot(v1);
					if(std::abs(1.0-d) > 1e-5 && std::abs(d) > 1e-5)
						break;
				}
				const auto c = v0.dot(v1);
				// 適当に倍率を掛ける
				v0 *= Scale<value_t>(mt);
				v1 *= Scale<value_t>(mt);

				struct {
					vec_t	v0, v1;
					value_t	cos;
				} ret{v0, v1, c};
				return ret;
			}
		}
		TYPED_TEST(VectorD_3, Cross) {
			using value_t = typename TestFixture::value_t;
			using vec_t = typename TestFixture::vec_t;

			const auto v = GenVecPair<vec_t>(this->mt());
			const auto s = std::sqrt(1 - v.cos*v.cos);
			const auto v2 = v.v0.cross(v.v1);
			// 外積の長さはv0とv1のそれを掛けた値*sinθ
			constexpr auto Th = Threshold<value_t>(0.7, 0);
			EXPECT_NEAR(v.v0.length()*v.v1.length()*s, v2.length(), Th);
			// 生成元のベクトルと内積を取ると0
			EXPECT_NEAR(v.v0.dot(v2), 0, Th);
			EXPECT_NEAR(v.v1.dot(v2), 0, Th);
			// operator % はcross()と同じ
			const auto v3 = v.v0 % v.v1;
			EXPECT_LE(AbsMax(vec_t(v2-v3)), Threshold<value_t>(0.3, 0));
		}
		TYPED_TEST(VectorD_3, VerticalVector) {
			using value_t = typename TestFixture::value_t;
			using vec_t = typename TestFixture::vec_t;
			const auto mtf = this->mt().template getUniformF<value_t>();
			const auto v0 = random::GenVecUnit<vec_t>(mtf);
			const auto v1 = v0.verticalVector();
			// ゼロベクトルではない
			EXPECT_GT(v1.length(), Threshold<value_t>(0.8,0));
			// 元のベクトルとの内積は0
			constexpr auto Th = Threshold<value_t>(0.2, 0);
			EXPECT_LE(std::abs(v0.dot(v1)), Th);
			// 生成したベクトルと外積を取ると、先の2つのベクトルとどちらも内積が0
			const auto v2 = (v1.cross(v0)).normalization();
			EXPECT_LE(std::abs(v2.dot(v0)), Th);
			EXPECT_LE(std::abs(v2.dot(v1)), Th);
		}
		// }
	}
}
