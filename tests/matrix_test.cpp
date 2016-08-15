#include "test.hpp"
#include "../compiler_macro.hpp"

namespace frea {
	namespace test {
		template <class T>
		using Matrix = RMatrix<T>;
		using Types_t = types::MatrixRange_t<types::Reg_t, 3,5, 3,5>;
		using Types = ToTestTypes_t<Types_t>;
		TYPED_TEST_CASE(Matrix, Types);

		// 内部表現による演算結果の差異をチェック
		TYPED_TEST(Matrix, Register) {
			using value_t = typename TestFixture::value_t;
			using array_t = typename TestFixture::array_t;

			// 逆数を掛けることで除算としていたりx86の場合内部精度がvalue_tと違うかも知れないのである程度のマージンを設ける
			constexpr auto Th = Threshold<value_t>(0.1, 2);		// 除算以外
			constexpr auto ThD = Threshold<value_t>(0.8, 2);	// 除算
			constexpr auto range = Range<value_t>{-1e3, 1e3};
			auto mat = this->makeRMat(range);
			const array_t raw(mat);
			const auto t = this->makeRMat(range);
			value_t s;
			auto rd = this->mt().template getUniformF<value_t>(range);
			do {
				s = rd();
			} while(IsZero(s, Th));
			#define DEF_TEST(num, op, target, th) \
				case num: { \
					const decltype(raw) raw0 = raw op target; \
					const decltype(mat) mat0 = mat op target; \
					/* 結果がほぼ一致することを確認 */ \
					ASSERT_LT(AbsMax(raw0 - mat0), th); \
					break; }

			// ランダムで四則演算
			switch(this->mt().template getUniform<int>({0,6})) {
				DEF_TEST(0, +, s, Th)		// スカラとの和
				DEF_TEST(1, -, s, Th)		// スカラとの差
				DEF_TEST(2, *, s, Th)		// スカラとの積
				DEF_TEST(3, /, s, ThD)		// スカラとの除
				DEF_TEST(4, +, t, Th)		// 行列との和
				DEF_TEST(5, -, t, Th)		// 行列との差
				DEF_TEST(6, *, t, Th)		// 行列との積
				default:
					__assume(false)
			}
			#undef DEF_TEST
		}
	}
}
