#include "test.hpp"

namespace frea {
	namespace test {
		// 整数ベクトルテストケース
		template <class T>
		using IntVector = RVector<T>;
		using ITypes = ToTestTypes_t<types::VectorRange_t<types::IReg_t>>;
		TYPED_TEST_CASE(IntVector, ITypes);

		// 整数の積算、除算に関するチェック
		TYPED_TEST(IntVector, MulDiv) {
			using value_t = typename TestFixture::value_t;
			using Lim = std::numeric_limits<value_t>;
			const Range<value_t> range = {Lim::lowest()/512, Lim::max()/512};
			const auto v0 = this->makeRVec(range);
			const int n = this->mt().template getUniform<int>({1,8});
			auto mul = v0 * n;
			auto sum = v0;
			for(int i=1 ; i<n ; i++)
				sum += v0;
			ASSERT_EQ(sum, mul);

			using vec_t = typename TestFixture::vec_t;
			// 自身との除算は1
			ASSERT_EQ(v0/v0, vec_t(1));
			{
				// 乗算して同じ数で除算すれば元と同じになる
				constexpr auto BitWidth = sizeof(value_t)*8;
				constexpr auto Width = (1 << BitWidth/4) - 1;
				const vec_t v0 = this->makeRVec({-Width, Width});
				const vec_t v1 = this->makeRVecNZ(1, {-Width, Width});
				ASSERT_EQ(v0, v0*v1/v1);
			}
		}
	}
}
