#include "test.hpp"
#include "../ulps.hpp"
#include "lubee/src/compiler_macro.hpp"

namespace frea {
	namespace test {
		template <class T>
		struct ULPs : Random {};
		using UTypes = ToTestTypes_t<types::Float_t>;
		TYPED_TEST_SUITE(ULPs, UTypes);

		TYPED_TEST(ULPs, Value) {
			using value_t = TypeParam;
			auto& mt = this->mt();
			const auto fvalue = mt.template getUniform<value_t>({std::numeric_limits<value_t>::lowest()/2,
																std::numeric_limits<value_t>::max()/2});
			const auto ivalueC = ulps::AsIntegral_C(fvalue);
			const auto ivalue = ulps::AsIntegral(fvalue);
			EXPECT_EQ(ivalueC, ivalue);
		}
		TYPED_TEST(ULPs, Compare) {
			using value_t = TypeParam;
			auto& mt = this->mt();
			auto fnRand = [&mt]() {
				const value_t range0 = std::pow(value_t(2), std::numeric_limits<value_t>::max_exponent/2);
				return mt.template getUniform<value_t>({-range0, range0});
			};

			const auto f0 = fnRand(),
						f1 = fnRand();
			// ULPs単位でどの程度離れているか
			const auto ulps = ulps::Diff(f0, f1),
						ulps_1 = std::max(decltype(ulps)(0), ulps-1);
			if(ulps < 0)
				return;

			// Equal範囲チェック
			EXPECT_TRUE(ulps::Equal(f0, f1, ulps));
			if(ulps > 0)
				EXPECT_FALSE(ulps::Equal(f0, f1, ulps_1));

			// LessEqual範囲チェック
			EXPECT_TRUE(ulps::LessEqual(f0, f1, ulps));
			if(ulps > 0) {
				if(f0 < f1)
					EXPECT_TRUE(ulps::LessEqual(f0, f1, ulps_1));
				else
					EXPECT_FALSE(ulps::LessEqual(f0, f1, ulps_1));
			}

			// GreaterEqual範囲チェック
			EXPECT_TRUE(ulps::GreaterEqual(f0, f1, ulps));
			if(ulps > 0) {
				if(f0 > f1)
					EXPECT_TRUE(ulps::GreaterEqual(f0, f1, ulps_1));
				else
					EXPECT_FALSE(ulps::GreaterEqual(f0, f1, ulps_1));
			}
		}
	}
}
