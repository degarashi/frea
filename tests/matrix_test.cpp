#include "test.hpp"
#include "../compiler_macro.hpp"
#include "../ulps.hpp"

namespace frea {
	namespace test {
		template <class T>
		using Matrix = RMatrix<T>;
		using Types_t = types::MatrixRange_t<types::Reg_t, 3,5, 3,5>;
		using Types = ToTestTypes_t<Types_t>;
		TYPED_TEST_CASE(Matrix, Types);

		TYPED_TEST(Matrix, Iterator) {
			USING(value_t);
			constexpr auto range = Range<value_t>{-1e3, 1e3};
			// イテレータを介して値を加算した場合とインデックスを介した場合で比べる　
			auto m0 = this->makeRMat(range);
			auto m1 = m0;
			const auto value = this->mt().template getUniform<value_t>(range);
			constexpr int M = decltype(m0)::dim_m,
						N = decltype(m1)::dim_n;
			for(auto itr=m1.begin() ; itr!=m1.end() ; ++itr)
				*itr += value;
			for(int i=0 ; i<M ; i++) {
				for(int j=0 ; j<N ; j++) {
					m0[i][j] += value;
				}
			}
			for(int i=0 ; i<M ; i++) {
				for(int j=0 ; j<N ; j++) {
					ASSERT_FLOAT_EQ(m0[i][j], m1[i][j]);
				}
			}

			// const_iteratorで巡回して値を比べる
			auto itr0 = m0.cbegin();
			auto itr1 = m1.cbegin();
			for(;;) {
				bool b0 = itr0==m0.cend(),
					 b1 = itr1==m1.cend();
				ASSERT_EQ(b0, b1);
				if(b0)
					break;
				ASSERT_FLOAT_EQ(*itr0, *itr1);
				++itr0;
				++itr1;
			}
		}
		TYPED_TEST(Matrix, IsZero) {
			USING(value_t);
			USING(mat_t);
			// 全ての要素がThreshold未満ならTrue
			const auto mtf = this->mt().template getUniformF<value_t>();
			const auto Th = mtf({ulps::Increment(value_t(0)), 1e3});
			const auto range = Range<value_t>{ulps::Increment(-Th), ulps::Decrement(Th)};
			auto m = this->makeRMat(range);
			ASSERT_TRUE(m.isZero(Th));

			// どこか一つ値をゼロ以外に変更すればFalse
			const auto idx = this->mt().template getUniformF<int>();
			const int idxM = idx({0, mat_t::dim_m-1}),
					idxN = idx({0, mat_t::dim_n-1});
			m[idxM][idxN] += Th*3;
			ASSERT_FALSE(m.isZero(Th));

			// 変更した行のみisZeroRowがFalse
			for(int i=0 ; i<mat_t::dim_m ; i++) {
				ASSERT_EQ(i!=idxM, m.isZeroRow(i, Th));
			}
		}
		TYPED_TEST(Matrix, Compare) {
			USING(mat_t);
			USING(value_t);
			constexpr auto range = Range<value_t>{-1e3, 1e3};
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
		TYPED_TEST(Matrix, Translation) {
			USING(mat_t);
			USING(vec_t);
			USING(value_t);
			constexpr auto range = Range<value_t>{-1e3, 1e3};
			auto v0 = this->makeRVec(range),
				vd = this->makeRVec(range);
			v0[vec_t::size-1] = 1;
			vd[vec_t::size-1] = 1;
			const auto m = mat_t::Translation(vd);
			const vec_t v1 = v0 * m;
			v0 += vd;
			v0[vec_t::size-1] = 1;
			EXPECT_LE(AbsMax(vec_t(v1-v0)), Threshold<value_t>(0.1, 0));
		}
		TYPED_TEST(Matrix, Scaling) {
			USING(mat_t);
			USING(vec_t);
			USING(value_t);
			using vecmin_t = typename vec_t::template type_cn<mat_t::dim_min>;
			using column_t = typename mat_t::column_t;

			constexpr auto range = Range<value_t>{-1e2, 1e2};
			const auto mtf = this->mt().template getUniformF<value_t>(range);
			const column_t v0 = random::GenVec<vecmin_t>(mtf).template convert<mat_t::dim_m>();
			const vecmin_t sc = random::GenVec<vecmin_t>(mtf);
			const auto m = mat_t::Scaling(sc);
			const vec_t v1a = v0 * m,
						v1b = v0 * sc;
			EXPECT_LE(AbsMax(vec_t(v1a-v1b)), Threshold<value_t>(0.1, 0));
		}
		// 内部表現による演算結果の差異をチェック
		TYPED_TEST(Matrix, Register) {
			USING(value_t);
			USING(array_t);

			// 逆数を掛けることで除算としていたりx86の場合内部精度がvalue_tと違うかも知れないのである程度のマージンを設ける
			constexpr auto Th = Threshold<value_t>(0.1, 2);		// 除算以外
			constexpr auto ThD = Threshold<value_t>(0.8, 2);	// 除算
			constexpr auto range = Range<value_t>{-1e3, 1e3};
			auto mat = this->makeRMat(range);
			const array_t raw(mat);
			const auto t = this->makeRMat(range);
			value_t s;
			auto rd = this->mt().template getUniformF<value_t>(range);
			do {
				s = rd();
			} while(IsZero(s, Th));
			#define DEF_TEST(num, op, target, th) \
				case num: { \
					const decltype(raw) raw0 = raw op target; \
					const decltype(mat) mat0 = mat op target; \
					/* 結果がほぼ一致することを確認 */ \
					ASSERT_LT(AbsMax(raw0 - mat0), th); \
					break; }

			// ランダムで四則演算
			switch(this->mt().template getUniform<int>({0,6})) {
				DEF_TEST(0, +, s, Th)		// スカラとの和
				DEF_TEST(1, -, s, Th)		// スカラとの差
				DEF_TEST(2, *, s, Th)		// スカラとの積
				DEF_TEST(3, /, s, ThD)		// スカラとの除
				DEF_TEST(4, +, t, Th)		// 行列との和
				DEF_TEST(5, -, t, Th)		// 行列との差
				DEF_TEST(6, *, t, Th)		// 行列との積
				default:
					__assume(false)
			}
			#undef DEF_TEST
		}
	}
}
