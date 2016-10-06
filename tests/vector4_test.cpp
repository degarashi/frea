#include "test.hpp"
#include "../lubee/compiler_macro.hpp"

namespace frea {
	namespace test {
		template <class T>
		using VectorD_4 = RVector<T>;
		using TypesD_4 = ToTestTypes_t<types::VectorRangeD_t<types::Float_t, 4>>;
		TYPED_TEST_CASE(VectorD_4, TypesD_4);

		namespace {
			template <class V, class RD, class T>
			V MakeWVec(RD&& rdf, const T& th) {
				V v;
				do {
					v = random::GenVec<V>(rdf);
				} while(v.w < th);
				return v;
			}
		}
		TYPED_TEST(VectorD_4, Wrap_Equality) {
			USING(vec_t);
			USING(value_t);
			{
				// asVec3Coord
				const auto v = MakeWVec<vec_t>(
					this->mt().template getUniformF<value_t>({-1e2, 1e2}),
					lubee::ThresholdF<value_t>(0.5)
				);
				const auto w = v.asInternal();
				EXPECT_EQ(v.asVec3Coord(), w.asVec3Coord());
			}
		}
		TYPED_TEST(VectorD_4, Vec3Coord) {
			USING(vec_t);
			USING(value_t);
			const auto v = MakeWVec<vec_t>(
				this->mt().template getUniformF<value_t>({-1e2, 1e2}),
				lubee::ThresholdF<value_t>(0.5)
			);
			using vec3_t = typename vec_t::template type_cn<3>;
			const vec3_t v0 = v.asVec3Coord();
			const vec3_t v1(v.x/v.w, v.y/v.w, v.z/v.w);
			EXPECT_LE(AbsMax(vec3_t(v0-v1)), lubee::ThresholdF<value_t>(0.8));
		}
	}
}
