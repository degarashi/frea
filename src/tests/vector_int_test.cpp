#include "test.hpp"

namespace frea {
	namespace test {
		// 整数ベクトルテストケース
		template <class T>
		using IntVector = RVector<T>;
		using ITypes = ToTestTypes_t<types::VectorRange_t<types::Int_t>>;
		TYPED_TEST_SUITE(IntVector, ITypes);

		// 整数の積算、除算に関するチェック
		TYPED_TEST(IntVector, MulDiv) {
			USING(value_t);
			USING(vec_t);
			{
				using Lim = std::numeric_limits<value_t>;
				const lubee::Range<value_t> range = {Lim::lowest()/2048, Lim::max()/2048};
				constexpr auto Th = lubee::Threshold<value_t>(0.4, 1);
				const auto v0 = this->makeRVecNZ(Th, range);
				const int n = this->mt().template getUniform<int>({1,8});
				auto mul = v0 * n;
				auto sum = v0;
				for(int i=1 ; i<n ; i++)
					sum += v0;
				ASSERT_EQ(sum, mul);
				// 自身との除算は1
				ASSERT_EQ(v0/v0, vec_t(1));
			}
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
