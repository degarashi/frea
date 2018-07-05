#include "test.hpp"
#include "../interpolation.hpp"

namespace frea {
	namespace test {
		// 浮動小数点数ベクトルテストケース
		template <class T>
		using FloatVector = RVector<T>;
		using FTypes = ToTestTypes_t<types::VectorRange_t<types::Float_t>>;
		TYPED_TEST_CASE(FloatVector, FTypes);

		namespace {
			template <class T>
			constexpr lubee::Range<T> DefaultRange{-1e3, 1e3};
		}
		TYPED_TEST(FloatVector, FunctionEquality_Method) {
			USING(vec_t);
			USING(value_t);
			const auto v = this->makeRVec(DefaultRange<value_t>);
			constexpr auto Th = lubee::ThresholdF<value_t>(0.1);
			{
				const auto w = v.asInternal();
				// Wrap
				{
					// normalize
					auto w0 = w.normalization(),
						 w1 = w;
					w1.normalize();
					EXPECT_LE(AbsMax(vec_t(w0-w1)), Th);
				}
				{
					// linearNormalize
					auto w0 = w.linearNormalization(),
						w1 = w;
					w1.linearNormalize();
					EXPECT_LE(AbsMax(vec_t(w0-w1)), Th);
				}
			}
			{
				// Vector
				{
					// normalize
					auto v0 = v.normalization(),
						 v1 = v;
					v1.normalize();
					EXPECT_LE(AbsMax(vec_t(v0-v1)), Th);
				}
				{
					// linearNormalize
					auto v0 = v.linearNormalization(),
						 v1 = v;
					v1.linearNormalize();
					EXPECT_LE(AbsMax(vec_t(v0-v1)), Th);
				}
			}
		}
		namespace {
			template <class V, class RD>
			V MakeNZVec(RD&& rdf, const typename V::value_t& th) {
				V ret;
				for(;;) {
					ret = random::GenVec<V>(rdf);
					bool b = true;
					for(int i=0 ; i<V::size ; i++) {
						if(std::abs(ret[i]) < th) {
							b = false;
							break;
						}
					}
					if(b)
						break;
				}
				return ret;
			}
		}
		TYPED_TEST(FloatVector, Wrap_Equality) {
			USING(value_t);
			USING(vec_t);
			const auto mtf = this->mt().template getUniformF<value_t>({-1e2, 1e2});
			auto v0 = MakeNZVec<vec_t>(mtf, 1e-2),
				v1 = MakeNZVec<vec_t>(mtf, 1e-2);
			auto w0 = v0.asInternal(),
				w1 = v1.asInternal();
			const auto mtf01 = this->mt().template getUniformF<value_t>({0,1});
			const auto s01 = mtf01();

			// distance
			EXPECT_EQ(v0.distance(v1), w0.distance(w1));
			// dist_sq
			EXPECT_EQ(v0.dist_sq(v1), w0.dist_sq(w1));
			// normalize
			EXPECT_EQ(vec_t(v0.normalization()), vec_t(w0.normalization()));
			// length
			EXPECT_EQ(v0.length(), w0.length());
			// len_sq
			EXPECT_EQ(v0.len_sq(), w0.len_sq());
			// l_intp
			EXPECT_EQ(v0.l_intp(v1, s01), w0.l_intp(w1, s01));
			// linearNormalization
			EXPECT_EQ(vec_t(v0.linearNormalization()), vec_t(w0.linearNormalization()));

			// isNaN & isOutstanding
			const auto chk = [&v0, &w0](){
				EXPECT_EQ(v0.isNaN(), w0.isNaN());
				EXPECT_EQ(v0.isOutstanding(), w0.isOutstanding());
			};
			using Lm = std::numeric_limits<value_t>;
			const int idx = this->mt().template getUniform<int>({0,vec_t::size-1});
			EXPECT_NO_FATAL_FAILURE(chk());
			v0[idx] = Lm::quiet_NaN();
			w0 = v0.asInternal();
			EXPECT_NO_FATAL_FAILURE(chk());
			v0[idx] = Lm::infinity();
			w0 = v0.asInternal();
			EXPECT_NO_FATAL_FAILURE(chk());
		}
		TYPED_TEST(FloatVector, LinearNormalize) {
			USING(value_t);
			USING(vec_t);
			const lubee::Range<value_t> range{-1e2, 1e2};
			const vec_t v = this->makeRVec(range),
						vn0 = v.linearNormalization();
			const value_t mv = v.absolute().getMaxValue();
			const vec_t vn1(vn0 * mv);
			// 割った数をかければ元と大体同じ
			EXPECT_LE(AbsMax(vec_t(vn1 - v)), lubee::ThresholdF<value_t>(0.3));
		}
		TYPED_TEST(FloatVector, Distance) {
			USING(vec_t);
			USING(array_t);
			USING(value_t);
			const lubee::Range<value_t> range{-1e2, 1e2};
			const auto v0 = this->makeRVec(range),
						v1 = this->makeRVec(range);
			const array_t a0(v0),
							a1(v1);
			value_t sum = 0;
			for(int i=0 ; i<vec_t::size ; i++)
				sum += Square(std::abs(a0[i] - a1[i]));
			constexpr auto Th = lubee::ThresholdF<value_t>(0.8);
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

			constexpr auto Th = lubee::ThresholdF<value_t>(0.1);
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
			// v2 = Vec::l_intp()にて計算
			// v3 = Vec::Internal::l_intp()にて計算
			const vec_t v2 = v0.l_intp(v1, t),
						v3 = v0.asInternal().l_intp(v1, t);
			// a2 = Arrayにて計算
			const array_t a2 = a0 + (a1 - a0) * t;
			// l2 = Lerp(Iterator)にて計算
			array_t l2;
			Lerp(l2.m, reinterpret_cast<const value_t*>(v0.m),
					reinterpret_cast<const value_t*>(v0.m + vec_t::size),
					reinterpret_cast<const value_t*>(v1.m), t);

			constexpr auto Th = lubee::ThresholdF<value_t>(0.7);
			// v2とv3の計算結果は同じ
			ASSERT_LT(AbsMax(vec_t(v3 - v2)), Th);
			// v2とa2も同じ
			ASSERT_LT(AbsMax(a2 - v2), Th);
			// l2とv2も同じ
			ASSERT_LT(AbsMax(l2 - v2), Th);
		}
		TYPED_TEST(FloatVector, MulDiv) {
			USING(value_t);
			constexpr auto threshold = lubee::ThresholdF<value_t>(0.7);
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
