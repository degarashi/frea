#include "test.hpp"

namespace frea {
	namespace test {
		// isNaN, Outstanding

		// 浮動小数点数ベクトルテストケース
		template <class T>
		using FloatVector = RVector<T>;
		using FTypes = ToTestTypes_t<FTypes_t>;
		TYPED_TEST_CASE(FloatVector, FTypes);

		TYPED_TEST(FloatVector, MulDiv) {
			using value_t = typename TestFixture::value_t;
			constexpr auto threshold = ThresholdULPs<value_t>;
			constexpr value_t rangev = RangeV<value_t>;
			const Range<value_t> range = {-rangev, rangev};
			const auto v0 = this->makeRVec(range);
			const int n = this->mt().template getUniform<int>({1,8});
			auto mul = v0 * n;
			auto sum = v0;
			for(int i=1 ; i<n ; i++)
				sum += v0;

			using vec_t = typename TestFixture::vec_t;
			{
				const vec_t tmul(mul),
							 tsum(sum);
				for(int i=0 ; i<TestFixture::size ; i++)
					ASSERT_NEAR(tsum.m[i], tmul.m[i], threshold);
			}
			// 自身との除算は1
			ASSERT_EQ(v0/v0, vec_t(1));
			// 乗算して同じ数で除算すれば大体元と同じになる
			const vec_t v1 = this->makeRVecNZ(1e-4, range);
			ASSERT_LT(DiffSum(v0.m, vec_t(v0*v1/v1).m), threshold);
		}
	}
}
