#include "test.hpp"
#include "../plane.hpp"
#include "../random/plane.hpp"

namespace frea {
	namespace test {
		template <class T>
		using VectorD_3 = RVector<T>;
		using TypesD_3 = ToTestTypes_t<types::VectorRangeD_t<types::FReg_t, 3>>;
		TYPED_TEST_CASE(VectorD_3, TypesD_3);

		TYPED_TEST(VectorD_3, Wrap_Equality) {
			USING(value_t);
			USING(vec_t);
			using quat_t = QuatT<value_t, false>;
			const auto mtf = this->mt().template getUniformF<value_t>();
			const auto dir = MakeDir<2, vec_t>(mtf, 0.8);
			const auto w0 = dir[0].asInternal(),
						w1 = dir[1].asInternal();
			const auto q = random::GenQuat<quat_t>(mtf);
			// cross
			EXPECT_EQ(dir[0].cross(dir[1]), w0.cross(w1));
			// op %
			EXPECT_EQ(dir[0]%dir[1], w0%w1);
			// verticalVector
			EXPECT_EQ(dir[0].verticalVector(), w0.verticalVector());
			// op * q
			EXPECT_EQ(dir[0]*q, w0*q);
		}
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
			USING(value_t);
			USING(vec_t);

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
			USING(value_t);
			USING(vec_t);
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
		TYPED_TEST(VectorD_3, Flip) {
			USING(value_t);
			USING(vec_t);
			constexpr Range<value_t> range{-1e2, 1e2};
			const auto mtf = this->mt().template getUniformF<value_t>(range);
			const auto p = random::GenPlane<PlaneT<value_t, false>>(mtf);
			const auto v0 = random::GenVec<vec_t>(mtf);

			// 一回反転すればplaneとの内積が符号反転
			const auto v1 = p.flip(v0);
			const auto d0 = p.dot(v0),
						d1 = p.dot(v1);
			constexpr auto Th = Threshold<value_t>(0.6, 0);
			EXPECT_NEAR(d0, -d1, Th);

			// 二回反転すると元に戻る
			const auto v2 = p.flip(v1);
			const auto d2 = p.dot(v2);
			EXPECT_NEAR(d0, d2, Th);

			// 反転前と後を結ぶ方向は平面の法線と同一、または正対する
			const vec_t dir = (v1 - v0).normalization();
			const auto dd = dir.dot(p.getNormal());
			EXPECT_LE(std::abs(dd)-1, Threshold<value_t>(0.8, 0));
		}
		TYPED_TEST(VectorD_3, CrossPoint) {
			USING(value_t);
			USING(vec_t);
			constexpr Range<value_t> range{-1e2, 1e2};
			const auto mtf = this->mt().template getUniformF<value_t>(range);
			const auto p = random::GenPlane<PlaneT<value_t, false>>(mtf);
			constexpr auto Th = Threshold<value_t>(0.6, 0);
			// 平面の両方側に位置するランダム頂点を生成
			vec_t vf, vb;
			// 表側に位置する頂点
			do {
				vf = random::GenVec<vec_t>(mtf);
			} while(p.dot(vf) <= Th);
			// 裏側に位置する頂点
			do {
				vb = random::GenVec<vec_t>(mtf);
			} while(p.dot(vb) >= -Th);
			const auto res0 = p.crosspoint(vf, vb);
			// 交差フラグはTrueである
			ASSERT_TRUE(res0.cross);
			// 交点は面上にある
			ASSERT_LE(std::abs(p.dot(res0.point)), Th);
			// 逆順で計算しても結果は同じ
			const auto res1 = p.crosspoint(vb, vf);
			ASSERT_TRUE(res1.cross);
			// 最初の交点と同じ位置
			ASSERT_LE(AbsMax(vec_t(res1.point - res0.point)), Th);
		}
	}
}
