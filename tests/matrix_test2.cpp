#include "test.hpp"
#include "../serialization/matrix.hpp"

namespace frea {
	namespace test {
		TYPED_TEST(Matrix, Compare) {
			USING(mat_t);
			USING(value_t);
			constexpr auto range = lubee::Range<value_t>{-1e3, 1e3};
			auto m0 = this->makeRMat(range),
				m1 = m0;
			EXPECT_EQ(m0, m0);
			EXPECT_EQ(m0, m1);
			const auto idx = this->mt().template getUniformF<int>();
			const int idxM = idx({0, mat_t::dim_m-1}),
					idxN = idx({0, mat_t::dim_n-1});
			m1[idxM][idxN] += 1;
			EXPECT_NE(m0, m1);
		}
		TYPED_TEST(Matrix, Scaling) {
			USING(mat_t);
			USING(vec_t);
			USING(value_t);
			using vecmin_t = typename vec_t::template type_cn<mat_t::dim_min>;
			using column_t = typename mat_t::column_t;
		
			constexpr auto range = lubee::Range<value_t>{-1e2, 1e2};
			const auto mtf = this->mt().template getUniformF<value_t>(range);
			const column_t v0 = random::GenVec<vecmin_t>(mtf).template convert<mat_t::dim_m>(),
							sc = random::GenVec<vecmin_t>(mtf).template convert<mat_t::dim_m>();
			const auto m = mat_t::Scaling(sc.template convert<mat_t::dim_min>());
			using vec2_t = typename vec_t::template type_cn<mat_t::dim_n>;
			const vec2_t v1a((v0 * m).template convert<mat_t::dim_n>()),
						v1b((v0 * sc).template convert<mat_t::dim_n>());
			EXPECT_LE(AbsMax(vec_t(v1a-v1b)), lubee::Threshold<value_t>(0.1, 0));
		}
		namespace {
			template <class TestFixture>
			void TranslationTest(TestFixture& self, std::true_type) {
				USING(mat_t);
				USING(vec_t);
				USING(value_t);
				constexpr auto range = lubee::Range<value_t>{-1e3, 1e3};
				const auto mtf = self.mt().template getUniformF<value_t>(range);
				constexpr auto dim_m = mat_t::dim_m;
				using column_t = typename mat_t::column_t;
				auto v0 = random::GenVec<column_t>(mtf),
					vd = random::GenVec<column_t>(mtf);
				v0[dim_m-1] = 1;
				vd[dim_m-1] = 1;
				const mat_t m = mat_t::Translation(vd.template convert<mat_t::dim_min>());
				const vec_t v1 = v0 * m;
				v0 += vd;
				v0[dim_m-1] = 1;
				auto v2 = v1.template convert<dim_m>();
				v2[dim_m-1] = 1;
				EXPECT_LE(AbsMax(column_t(v2-v0)), lubee::Threshold<value_t>(0.1, 0));
			}
			template <class TestFixture>
			void TranslationTest(TestFixture&, std::false_type) {}
		}
		TYPED_TEST(Matrix, Translation) {
			USING(mat_t);
			ASSERT_NO_FATAL_FAILURE(TranslationTest(*this, lubee::BConst<mat_t::dim_m == mat_t::dim_n+1>()));
		}
		TYPED_TEST(Matrix, Serialization) {
			lubee::CheckSerialization(this->makeRMat());
		}
	}
}
