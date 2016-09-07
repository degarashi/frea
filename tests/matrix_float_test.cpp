#include "test.hpp"

namespace frea {
	namespace test {
		template <class T>
		using FMatrix = RMatrix<T>;
		using Types_t = types::MatrixRange_t<types::Float_t, 3,5, 3,5>;
		using Types = ToTestTypes_t<Types_t>;
		TYPED_TEST_CASE(FMatrix, Types);

		TYPED_TEST(FMatrix, LinearNormalize) {
			USING(value_t);
			USING(mat_t);
			USING(vec_t);

			constexpr auto range = Range<value_t>{-1e2, 1e2};
			const auto m = this->makeRMat(range),
						mn0 = m.linearNormalization();
			value_t r[mat_t::dim_m];
			for(int i=0 ; i<mat_t::dim_m ; i++)
				r[i] = m.m[i].absolute().getMaxValue();
			// 割った数をかければ元と大体同じ(行ごとに比較する)
			for(int i=0 ; i<mat_t::dim_m ; i++) {
				const vec_t& v0 = m.m[i];
				const vec_t v1 = mn0.m[i] * r[i];
				EXPECT_LE(AbsMax(vec_t(v1-v0)), ThresholdF<value_t>(0.3));
			}
		}
		TYPED_TEST(FMatrix, RowReduce) {
			USING(mat_t);
			USING(value_t);

			constexpr auto Th = ThresholdF<value_t>(0.1);
			constexpr auto range = Range<value_t>{1e2};
			constexpr auto dim_m = mat_t::dim_m,
						dim_n = mat_t::dim_n;
			{
				// RowReduceしたら被約階段行列になっている
				auto m = this->makeRMat(range);
				const int zr = m.rowReduce(Th);
				ASSERT_TRUE(m.isEchelon(Th));
				int count = 0;
				for(int i=0 ; i<dim_m ; i++)
					count += m[i].isZero(Th);
				// 戻り値とゼロ行の数は一致する
				ASSERT_EQ(zr, count);
			}

			if(dim_n == dim_m+1) {
				// 連立一次方程式を解いてみる
				auto m = this->makeRMat(range);
				const auto mtf = this->mt().template getUniformF<value_t>(Range<value_t>{1e1});
				// 変数値を適当に決める
				value_t var[dim_m];
				for(int i=0 ; i<dim_m ; i++)
					var[i] = mtf();
				// 変数値に合わせた計算結果を代入
				for(int i=0 ; i<dim_m ; i++) {
					value_t sum = 0;
					for(int j=0 ; j<dim_n-1 ; j++)
						sum += m[i][j] * var[j];
					m[i][dim_n-1] = sum;
				}

				const int zr = m.rowReduce(Th);
				if(zr == 0) {
					constexpr auto Th = ThresholdF<value_t>(0.9);
					// 解が一意に定まる場合はチェックする
					for(int i=0 ; i<dim_m ; i++)
						EXPECT_NEAR(m[i][dim_n-1], var[i], Th);
				}
			}
		}
	}
}
