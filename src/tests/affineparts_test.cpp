#include "test.hpp"
#include "../affine_parts.hpp"
#include "../serialization/affine_parts.hpp"

namespace frea {
	namespace test {
		template <class T>
		struct AffineParts : Random {
			using value_t = T;
			using AP = ::frea::AffineParts<value_t>;

			auto makeAP(const lubee::Range<value_t>& rt, const lubee::Range<value_t>& rs) {
				AP ap;
				using Vec = Vec_t<value_t, 3, true>;
				ap.offset = random::GenVec<Vec>(mt().template getUniformF<value_t>(rt));
				ap.scale = random::GenVec<Vec>(mt().template getUniformF<value_t>(rs));
				using Quat = QuatT<value_t, true>;
				for(;;) {
					ap.rotation = random::GenQuat<Quat>(mt().template getUniformF<value_t>());
					try {
						ap.rotation.getAxis();
						break;
					} catch(const NoValidAxis&) {}
				}
				return ap;
			}
			auto makeAP() {
				const lubee::Range<value_t> r{-1e3, 1e3};
				return makeAP(r,r);
			}
		};
		using Types = ToTestTypes_t<types::Float_t>;
		TYPED_TEST_SUITE(AffineParts, Types);

		TYPED_TEST(AffineParts, Test) {
			// 行列をDecompAffineした結果を再度合成して同じかどうかチェック
			const auto ap = this->makeAP({-1e3, 1e3}, {1e-2, 1e2});

			USING(value_t);
			using M33 = Mat_t<value_t, 3,3, true>;
			using M44 = Mat_t<value_t, 4,4, true>;
			const auto mR = M33::RotationAxis(ap.rotation.getAxis(), ap.rotation.angle()).template convertI<4,4>(1);
			const auto mS = M44::Scaling(ap.scale.template convertI<4,3>(1));
			const auto mT = M44::Translation(ap.offset.template convertI<4,3>(1));
			const auto m = mS * mR * mT;

			USING(AP);
			const auto ap2 = AP::Decomp(m);
			using V3 = Vec_t<value_t, 3, true>;

			const auto mtf = this->mt().template getUniformF<value_t>();
			auto dir0 = random::GenVecUnit<V3>(mtf),
				dir1 = dir0;
			dir0 *= ap.rotation;
			dir1 *= ap2.rotation;
			constexpr auto Th = lubee::ThresholdF<value_t>(0.5);
			EXPECT_LT(AbsMax(V3(ap.offset - ap2.offset)), Th);
			EXPECT_LT(AbsMax(V3(ap.scale - ap2.scale)), Th);
			EXPECT_LT(AbsMax(V3(dir0 - dir1)), Th);
		}
		TYPED_TEST(AffineParts, Serialization) {
			lubee::CheckSerialization(this->makeAP());
		}
	}
}
