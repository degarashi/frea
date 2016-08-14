#include "test.hpp"
#include "../compiler_macro.hpp"

namespace frea {
	namespace test {
		// template <class T>
		// using SMatrix = RMatrix<T>;
		// using SqTypes_t = seq::ExpandTypes_t2<
		// 	types::SquareMat_t,
		// 	std::tuple<
		// 		types::Reg_t,
		// 		seq::Range_t<2,5>,
		// 		std::tuple<BConst<false>>
		// 	>
		// >;
		// using SqTypes = ToTestTypes_t<SqTypes_t>;
		// TYPED_TEST_CASE(SMatrix, SqTypes);
		// TYPED_TEST(SMatrix, Transpose) {
		// 	using value_t = typename TestFixture::value_t;
		// 	using array_t = typename TestFixture::array_t;
		// 	using mat_t = typename TestFixture::mat_t;
		//
		// 	// 転置行列を二次元配列で計算した結果と比較
		// 	constexpr auto range = Range<value_t>{-1e3, 1e3};
		// 	auto mat = this->makeRMat(range);
		// 	array_t raw(mat);
		//
		// 	raw.transpose();
		// 	mat.transpose();
		// 	ASSERT_LT(AbsMax(raw - mat), Threshold<value_t>());
		//
		// 	// (AB)t と (B)t(A)tの結果は同じ
		// 	mat = this->makeRMat(range);
		// 	auto mat2 = this->makeRMat(range);
		// 	const mat_t result0 = (mat * mat2).transposition();
		// 	mat.transpose();
		// 	mat2.transpose();
		// 	const mat_t result1 = mat2 * mat;
		// 	ASSERT_LT(AbsMax(mat_t(result0 - result1)), Threshold<value_t>());
		// }
	}
}
