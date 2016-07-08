#include "test.hpp"
#include "../compiler_macro.hpp"

namespace frea {
	namespace test {
		template <class R, class S, class A>
		using SquareMat_t = std::tuple<R, S, S, A>;
		using SqTypes_t = seq::ExpandTypes_t2<
			SquareMat_t,
			std::tuple<
				std::tuple<__m128, __m128d>,
				seq::Range_t<2,5>,
                seq::BoolSeq_t
			>
		>;
		using SqTypes = ToTestTypes_t<SqTypes_t>;

		template <class T>
		using FMatrix = RMatrix<T>;
		TYPED_TEST_CASE(FMatrix, SqTypes);

		// 内部表現による演算結果の差異をチェック
		TYPED_TEST(FMatrix, Register) {
			using value_t = typename TestFixture::value_t;
			constexpr auto range = Range<value_t>{-1e3, 1e3};
			auto mat = this->makeRMat(range);
			typename TestFixture::array_t raw(mat);
			const auto t = this->makeRMat(range);
			value_t s;
			do {
				s = this->rdf()(range);
			} while(std::abs(s) < 1e-4);
			#define DEF_TEST(num, op, target) \
				case num: { \
					const decltype(raw) raw0 = raw op target; \
					const decltype(mat) mat0 = mat op target; \
					/* 結果がほぼ一致することを確認 */ \
					ASSERT_LT(AbsMax(raw0 - mat0), 1e-3); \
					break; }

			// ランダムで四則演算
			switch(this->mt().template getUniform<int>({0,6})) {
				DEF_TEST(0, +, s)		// スカラとの和
				DEF_TEST(1, -, s)		// スカラとの差
				DEF_TEST(2, *, s)		// スカラとの積
				DEF_TEST(3, /, s)		// スカラとの除
				DEF_TEST(4, +, t)		// 行列との和
				DEF_TEST(5, -, t)		// 行列との差
				DEF_TEST(6, *, t)		// 行列との積
				default:
					__assume(false)
			}
		}
	}
}
