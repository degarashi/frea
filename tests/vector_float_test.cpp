#include "test.hpp"

namespace frea {
	namespace test {
		// 浮動小数点数ベクトルテストケース
		template <class T>
		using FloatVector = RVector<T>;
		using FTypes = ToTestTypes_t<types::VectorRange_t<types::FReg_t>>;
		TYPED_TEST_CASE(FloatVector, FTypes);

		namespace {
			template <class T>
			constexpr Range<T> DefaultRange{-1e3, 1e3};
		}
		TYPED_TEST(FloatVector, Distance) {
			USING(vec_t);
			USING(array_t);
			USING(value_t);
			const Range<value_t> range{-1e2, 1e2};
			const auto v0 = this->makeRVec(range),
						v1 = this->makeRVec(range);
			const array_t a0(v0),
							a1(v1);
			value_t sum = 0;
			for(int i=0 ; i<vec_t::size ; i++)
				sum += Square(std::abs(a0[i] - a1[i]));
			constexpr auto Th = ThresholdF<value_t>(0.8);
			ASSERT_NEAR(sum, v0.dist_sq(v1), Th);
			sum = std::sqrt(sum);
			ASSERT_NEAR(sum, v0.distance(v1), Th);
		}
		TYPED_TEST(FloatVector, Normalize) {
			USING(vec_t);
			USING(value_t);

			vec_t vec = this->makeRVec(DefaultRange<value_t>);
			auto iv = vec.asInternal();
			iv.normalize();
			vec.normalize();

			constexpr auto Th = ThresholdF<value_t>(0.1);
			ASSERT_LE(AbsMax(vec_t(iv-vec)), Th);
			ASSERT_NEAR(vec.length(), value_t(1.0), Th);
		}
		TYPED_TEST(FloatVector, Interpolation) {
			USING(vec_t);
			USING(value_t);
			USING(array_t);

			const vec_t v0 = this->makeRVec(DefaultRange<value_t>),
						v1 = this->makeRVec(DefaultRange<value_t>);
			const array_t a0(v0),
						a1(v1);
			const value_t t = this->mt().template getUniform<value_t>({0,1});
			const vec_t v2 = v0.l_intp(v1, t),
						v3 = v0.asInternal().l_intp(v1, t);
			const array_t a2 = a0 + (a1 - a0) * t;

			constexpr auto Th = ThresholdF<value_t>(0.7);
			ASSERT_LT(AbsMax(vec_t(v3 - v2)), Th);
			ASSERT_LT(AbsMax(a2 - v2), Th);
		}
		TYPED_TEST(FloatVector, MulDiv) {
			USING(value_t);
			constexpr auto threshold = ThresholdF<value_t>(0.7);
			const auto v0 = this->makeRVec(DefaultRange<value_t>);
			const int n = this->mt().template getUniform<int>({1,8});
			auto mul = v0 * n;
			auto sum = v0;
			for(int i=1 ; i<n ; i++)
				sum += v0;

			USING(vec_t);
			{
				const vec_t tmul(mul),
							 tsum(sum);
				for(int i=0 ; i<TestFixture::size ; i++)
					ASSERT_NEAR(tsum.m[i], tmul.m[i], threshold);
			}
			// 自身との除算は1
			ASSERT_EQ(vec_t(v0/v0), vec_t(1));
			// 乗算して同じ数で除算すれば大体元と同じになる
			const vec_t v1 = this->makeRVecNZ(1e-4, DefaultRange<value_t>);
			ASSERT_LT(DiffSum(v0.m, vec_t(v0*v1/v1).m), threshold);
		}

		// isNaN, Outstandingのチェック
		TYPED_TEST(FloatVector, InvalidValue) {
			USING(value_t);
			USING(vec_t);
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
