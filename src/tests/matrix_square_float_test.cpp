#include "test.hpp"

namespace frea {
	namespace test {
		template <class T>
		using FSMatrix = RMatrix<T>;
		using FSTypes_t = types::SMatrixRange_t<types::Float_t, 2,5>;
		using FSTypes = ToTestTypes_t<FSTypes_t>;
		TYPED_TEST_SUITE(FSMatrix, FSTypes);
		namespace {
			template <class M>
			bool IsNormalMatrix(const M& m, const typename M::value_t& th) {
				for(int i=0 ; i<M::dim_m ; i++) {
					for(int j=0 ; j<M::dim_n ; j++) {
						const typename M::value_t val = std::abs(m[i][j] - ((i==j) ? 1 : 0));
						if(val >= th)
							return false;
					}
				}
				return true;
			}
		}
		TYPED_TEST(FSMatrix, Inverse) {
			USING(value_t);
			USING(mat_t);
			constexpr auto Th = lubee::ThresholdF<value_t>(0.9);

			constexpr auto range = lubee::Range<value_t>{-1, 1};
			const auto mat = this->makeRMat(range);
			try {
				const auto mati = mat.inversion();
				const mat_t m0(mat*mati),
							m1(mati*mat);
				ASSERT_TRUE(IsNormalMatrix(m0, Th));
				ASSERT_TRUE(IsNormalMatrix(m1, Th));
			} catch(const NoInverseMatrix&) {
				// 逆行列が存在しない場合は何もチェックしない
			}
		}
	}
}
