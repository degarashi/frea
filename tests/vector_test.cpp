#include "test.hpp"
#include "../compiler_macro.hpp"

namespace frea {
	namespace test {
		// operator -
		// mask(LH), setAt+initAt+pickAt, makeEquality
		// dot, average, distance, normalize, length
		// saturation
		// interpolation
		// convert, convertI

		// Int, Float両方のテストケース
		template <class T>
		using Vector = RVector<T>;
		using Types = ToTestTypes_t<types::VectorRange_t<types::Reg_t>>;
		TYPED_TEST_CASE(Vector, Types);

		template <class T>
		constexpr Range<T> DefaultRange{-1e3, 1e3};

		TYPED_TEST(Vector, SumUp) {
			using value_t = typename TestFixture::value_t;
			using array_t = typename TestFixture::array_t;
			constexpr auto range = DefaultRange<value_t>;
			const auto vec = this->makeRVec(range);
			const array_t ar(vec);

			value_t sum0 = 0;
			for(int i=0 ; i<array_t::size ; i++)
				sum0 += ar[i];
			const value_t sum1 = (vec+0).sumUp();
			constexpr auto Th = ulps::Diff_C<value_t>(0, Threshold<value_t>(1<<5, 0));
			ASSERT_TRUE(ulps::Equal(sum0, sum1, Th));
		}
		// 内部表現による演算結果の差異をチェック
		TYPED_TEST(Vector, Register) {
			using value_t = typename TestFixture::value_t;
			auto& mt = this->mt();
			constexpr auto range = DefaultRange<value_t>;
			auto vec = this->makeRVec(range);
			constexpr int size = decltype(vec)::size;
			typename TestFixture::array_t	raw;
			for(int i=0 ; i<size ; i++)
				raw.m[i] = vec.m[i];
			const auto chk = [&](auto dist){
				for(int i=0 ; i<size ; i++) {
					if(TestFixture::integral)
						ASSERT_NEAR(raw.m[i], vec.m[i], dist);
					else
						ASSERT_FLOAT_EQ(raw.m[i], vec.m[i]);
				}
			};
			const auto rdi = mt.template getUniformF<int>();
			const auto t = this->makeRVecNZ(value_t(TestFixture::integral ? 1 : 1e-4), {-1e4, 1e4});
			// ランダムで四則演算
			switch(rdi({0,3})) {
				case 0:
					vec += t;
					raw += t;
					// 結果が一致することを確認
					ASSERT_NO_FATAL_FAILURE(chk(0));
					break;
				case 1:
					vec -= t;
					raw -= t;
					ASSERT_NO_FATAL_FAILURE(chk(0));
					break;
				case 2:
					vec *= t;
					raw *= t;
					ASSERT_NO_FATAL_FAILURE(chk(0));
					break;
				case 3:
					vec /= t;
					raw /= t;
					ASSERT_NO_FATAL_FAILURE(chk(0));
					break;
				default:
					__assume(false)
			}
		}
		// 要素比較チェック
		TYPED_TEST(Vector, Compare) {
			using value_t = typename TestFixture::value_t;
			constexpr auto range = DefaultRange<value_t>;
			const auto v0 = this->makeRVec(range),
						v1 = this->makeRVec(range);
			// 自身と比較したら==はtrueになる
			ASSERT_EQ(v0, v0);
			// ==と!=は正反対の結果になる
			ASSERT_NE(v0==v1, v0!=v1);

			using vec_t = typename TestFixture::vec_t;
			const auto randNumNZ = [this]{
				auto ret = this->mt().template getUniform<value_t>();
				if(std::abs(ret) < 1)
					ret += 10;
				return ret;
			};
			// 自身を足し合わせれば値が変わる
			vec_t tv0 = v0 + v0;
			ASSERT_NE(tv0, v0);
			for(int i=0 ; i<vec_t::size ; i++) {
				tv0 = v0;
				ASSERT_EQ(tv0, v0);
				// 要素を1つだけいじっても==がfalseになる
				tv0.m[i] += randNumNZ();
				ASSERT_NE(tv0, v0);
			}
		}
		// 論理演算チェック
		TYPED_TEST(Vector, Logical) {
			using value_t = typename TestFixture::value_t;
			using vec_t = typename TestFixture::vec_t;
			constexpr auto range = DefaultRange<value_t>;
			const auto v0 = this->makeRVec(range);
			ASSERT_EQ(v0.asInternal(), v0);
			// 自身とのAndは元と同じ
			ASSERT_EQ(v0&v0, v0);
			// 自身とのOrも同じ
			ASSERT_EQ(v0|v0, v0);
			// 自身とのXorは必ずゼロになる
			ASSERT_EQ(vec_t(0), v0^v0);

			// (整数値のみ)
			if(TestFixture::integral) {
				const vec_t one = vec_t::I::One();
				// ~0とのAndは元と同じ
				ASSERT_EQ(v0, v0&one);
				// ~0とのOrは~0になる
				ASSERT_EQ(one, v0|one);
				// ~0のXorで2回反転させれば元と同じになる
				ASSERT_EQ(v0, (v0^one^one));
			}
		}
		// 四則演算のチェック
		TYPED_TEST(Vector, Arithmetic) {
			using value_t = typename TestFixture::value_t;
			constexpr auto range = DefaultRange<value_t>;
			const auto v0 = this->makeRVec(range),
						v1 = this->makeRVec(range);
			const auto check = [](auto v0, auto v1, const auto& op, const auto& ope) {
				// 代入後の値チェック
				auto tv0 = v0;
				const auto tv1 = v1;
				ASSERT_EQ(v0, tv0);
				ASSERT_EQ(v1, tv1);

				// 三項演算をどこにも代入しなければ値は変わらない
				op(v0, v1);
				ASSERT_EQ(tv0, v0);
				ASSERT_EQ(tv1, v1);
				// 二項演算と三項演算は一緒
				tv0 = op(v0, v1);
				ope(v0, v1);
				ASSERT_EQ(tv0, v0);
				// 複数の三項演算を二項演算に分けても結果は一緒
				v0 = tv0;
				v0 = op(op(v0, v0), v1);
				ope(tv0, tv0);
				ope(tv0, v1);
				ASSERT_EQ(tv0, v0);
			};
			#define DEF_TEST(op) \
				ASSERT_NO_FATAL_FAILURE( \
					check( \
						v0, v1, \
						[](const auto& p0, const auto& p1){ return p0 op p1; }, \
						[](auto& p0, const auto& p1){ return p0 op##= p1; } \
					) \
				);
			DEF_TEST(+)
			DEF_TEST(-)
			DEF_TEST(*)
			// Divideの場合は0除算を避ける
			if(!HasZero(v1.m, (TestFixture::integral ? 1: 1e-3f))) {
				DEF_TEST(/)
			}
			#undef DEF_TEST
		}
	}
}
