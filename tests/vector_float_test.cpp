#include "test.hpp"

namespace frea {
	namespace test {
		// 浮動小数点数ベクトルテストケース
		template <class T>
		using FloatVector = RVector<T>;
		using FTypes = ToTestTypes_t<types::VectorRange_t<types::FReg_t>>;
		TYPED_TEST_CASE(FloatVector, FTypes);

		TYPED_TEST(FloatVector, MulDiv) {
			using value_t = typename TestFixture::value_t;
			constexpr auto threshold = Threshold<value_t>(1<<14, 0);
			constexpr value_t rangev = 1e3;
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

		// isNaN, Outstandingのチェック
		TYPED_TEST(FloatVector, InvalidValue) {
			using value_t = typename TestFixture::value_t;
			using vec_t = typename TestFixture::vec_t;
			using VL = std::numeric_limits<value_t>;
			vec_t qn = this->makeRVec(),
				sn = this->makeRVec(),
				inf = this->makeRVec();
			auto& mt = this->mt();
			ASSERT_FALSE(qn.isNaN());
			ASSERT_FALSE(qn.isOutstanding());
			ASSERT_FALSE(sn.isNaN());
			ASSERT_FALSE(sn.isOutstanding());
			ASSERT_FALSE(inf.isNaN());
			ASSERT_FALSE(inf.isOutstanding());

			const int idx = mt.template getUniform<int>({0,vec_t::size-1});
			qn[idx] = VL::quiet_NaN();
			sn[idx] = VL::signaling_NaN();
			inf[idx] = VL::infinity();

			ASSERT_TRUE(qn.isNaN());
			ASSERT_TRUE(qn.isOutstanding());
			ASSERT_TRUE(sn.isNaN());
			ASSERT_TRUE(sn.isOutstanding());
			ASSERT_FALSE(inf.isNaN());
			ASSERT_TRUE(inf.isOutstanding());
		}
	}
}
