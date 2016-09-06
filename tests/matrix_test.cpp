#include "test.hpp"
#include "../compiler_macro.hpp"
#include "../ulps.hpp"

namespace frea {
	namespace test {
		namespace check {
			template <class M, class W>
			void GetRow(const M&, const W&, IConst<-1>) {}
			template <class M, class W, int N>
			void GetRow(const M& m, const W& w, IConst<N>) {
				ASSERT_EQ(m.template getRow<N>(), w.template getRow<N>());
				ASSERT_NO_FATAL_FAILURE(GetRow(m, w, IConst<N-1>()));
			}
			template <class M, class W>
			void GetColumn(const M&, const W&, IConst<-1>) {}
			template <class M, class W, int N>
			void GetColumn(const M& m, const W& w, IConst<N>) {
				ASSERT_EQ(m.template getColumn<N>(), w.template getColumn<N>());
				ASSERT_NO_FATAL_FAILURE(GetColumn(m, w, IConst<N-1>()));
			}

			template <class M, class W, class V>
			void SetRow(M&, W&, const V&, IConst<-1>) {}
			template <class M, class W, class V, int N>
			void SetRow(M& m, W& w, const V& r, IConst<N>) {
				const auto tmp0 = m.template getRow<N>();
				const auto tmp1 = w.template getRow<N>();
				m.template setRow<N>(r);
				w.template setRow<N>(r);
				ASSERT_EQ(m, w);
				m.template setRow<N>(tmp0);
				w.template setRow<N>(tmp1);
				ASSERT_NO_FATAL_FAILURE(SetRow(m, w, r, IConst<N-1>()));
			}
			template <class M, class W, class V>
			void SetColumn(M&, W&, const V&, IConst<-1>) {}
			template <class M, class W, class V, int N>
			void SetColumn(M& m, W& w, const V& c, IConst<N>) {
				const auto tmp0 = m.template getColumn<N>();
				const auto tmp1 = w.template getColumn<N>();
				m.template setColumn<N>(c);
				w.template setColumn<N>(c);
				ASSERT_EQ(m, w);
				m.template setColumn<N>(tmp0);
				w.template setColumn<N>(tmp1);
				ASSERT_NO_FATAL_FAILURE(SetColumn(m, w, c, IConst<N-1>()));
			}
		}
		TYPED_TEST(Matrix, Wrap_Equality) {
			USING(value_t);
			USING(mat_t);
			USING(vec_t);
			using wrap_t = typename mat_t::wrap_t;
			using column_t = typename mat_t::column_t;
			constexpr auto range = Range<value_t>{1e3};
			mat_t m0 = this->makeRMat(range),
					m1 = this->makeRMat(range);
			const auto mtf = this->mt().template getUniformF<value_t>(range);
			using mat2_t = typename mat_t::template type_cn<mat_t::dim_n, mat_t::dim_n>;
			const mat2_t m2 = random::GenMat<mat2_t>(mtf);
			auto w0 = m0.asInternal(),
					w1 = m1.asInternal();
			const auto w2 = m2.asInternal();
			const auto s = mtf();
			const column_t vc = random::GenVec<column_t>(mtf);
			const vec_t vr = random::GenVec<vec_t>(mtf);

			// (mat +-* mat) == (wrapM +-* wrapM)
			EXPECT_EQ(m0+m1, w0+w1);
			EXPECT_EQ(m0-m1, w0-w1);
			EXPECT_EQ(m0*m2, w0*w2);
			// (mat +-*/ s) == (wrapM +-*/ s)
			EXPECT_EQ(m0+s, w0+s);
			EXPECT_EQ(m0-s, w0-s);
			EXPECT_EQ(m0*s, w0*s);
			EXPECT_EQ(m0/s, w0/s);
			// (mat == mat) == (wrapM == wrapM)
			EXPECT_EQ(m0==m0, w0==w0);
			EXPECT_EQ(m0!=m0, w0!=w0);
			// Identity
			EXPECT_EQ(mat_t::Identity(), wrap_t::Identity());
			// (mat * vec) == (wrapM * vec)
			EXPECT_EQ(m0*vr, w0*vr);
			// (vec * mat) == (vec * wrapM)
			EXPECT_EQ(vc*m0, vc*w0);
			// diagonal
			EXPECT_EQ(mat_t::Diagonal(s), wrap_t::Diagonal(s));
			// linearNormalization
			EXPECT_EQ(m0.linearNormalization(), w0.linearNormalization());

			check::GetRow(m0, w0, IConst<mat_t::dim_m-1>());
			check::GetColumn(m0, w0, IConst<mat_t::dim_n-1>());
			check::SetRow(m0, w0, vr, IConst<mat_t::dim_m-1>());
			check::SetColumn(m0, w0, vc, IConst<mat_t::dim_n-1>());
		}
		TYPED_TEST(Matrix, FunctionEquality_Method) {
			USING(value_t);
			constexpr auto range = Range<value_t>{1e3};
			// linearNormalize
			auto m0 = this->makeRMat(range),
				 m1 = m0.linearNormalization();
			m0.linearNormalize();
			EXPECT_EQ(m0, m1);
		}
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
		TYPED_TEST(Matrix, IsEchelon) {
			USING(mat_t);
			USING(value_t);

			constexpr auto dim_m = mat_t::dim_m,
							dim_n = mat_t::dim_n;
			mat_t m;
			const auto rI = this->mt().template getUniformF<int>();
			constexpr auto range = Range<value_t>{1e2};
			const auto rV = this->mt().template getUniformF<value_t>(range);
			const auto th = Threshold<value_t>(0.1, 1);
			int edge[mat_t::dim_m];
			int curN = 0;
			for(int curM=0 ; curM<dim_m ; curM++) {
				edge[curM] = curN;
				// [curN-1]まではゼロで埋める
				for(int i=0 ; i<curN ; i++)
					m[curM][i] = 0;
				// [curM][curN] -> 1
				if(curN < dim_n)
					m[curM][curN] = 1;
				// [curN+1]から末尾までは適当に値を入れておく
				for(int i=curN+1 ; i<dim_n ; i++)
					m[curM][i] = rV();

				if(curN < dim_n) {
					// 横カーソルをランダムに進める
					const int adv = rI({1, dim_n-curN});
					curN += adv;
				}
			}
			ASSERT_TRUE(m.isEchelon(th));

			// 階段の縁より右側なら変更しても結果は変わらない
			for(int idxM=0 ; idxM<dim_m ; idxM++) {
				if(dim_n - edge[idxM] - 1 > 0) {
					const int idxN = rI({1, dim_n-edge[idxM]-1})+edge[idxM];
					m[idxM][idxN] = rV();
					EXPECT_TRUE(m.isEchelon(th));
				}
			}
			// 階段の縁を1以外にするとFalseになる
			for(int idxM=0 ; idxM<dim_m ; idxM++) {
				if(edge[idxM] < dim_n) {
					auto tm = m;
					value_t num;
					do { num=rV(); } while(std::abs(num-1) <= th);
					tm[idxM][edge[idxM]] = num;
					EXPECT_FALSE(tm.isEchelon(th));
				}
			}
			// 階段の縁から左を0以外にするとFalseになる
			for(int idxM=1 ; idxM<dim_m ; idxM++) {
				const int em = edge[idxM-1];
				if(em < dim_n) {
					const int idxN = rI({0, em});
					auto tm = m;
					value_t num;
					do { num=rV(); } while(std::abs(num) <= th);
					tm[idxM][idxN] = num;
					EXPECT_FALSE(tm.isEchelon(th));
				}
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
		// 内部表現による演算結果の差異をチェック
		TYPED_TEST(Matrix, Register) {
			USING(value_t);
			USING(array_t);
			USING(mat_t);
		
			// 逆数を掛けることで除算としていたりx86の場合内部精度がvalue_tと違うかも知れないのである程度のマージンを設ける
			constexpr auto Th = Threshold<value_t>(0.1, 2);		// 除算以外
			constexpr auto ThD = Threshold<value_t>(0.8, 2);	// 除算
			constexpr auto range = Range<value_t>{-1e3, 1e3};
			auto mat = this->makeRMat(range);
			const array_t raw(mat);

			const auto rd = this->mt().template getUniformF<value_t>(range);
			// 計算対象: 積算用
			using matM_t = typename mat_t::template type_cn<mat_t::dim_n, mat_t::dim_n>;
			const auto t_mul = random::GenMat<matM_t>(rd);
			// 計算対象: 加減算用
			const auto t_add = this->makeRMat(range);
			value_t s;
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
				DEF_TEST(4, +, t_add, Th)	// 行列との和
				DEF_TEST(5, -, t_add, Th)	// 行列との差
				DEF_TEST(6, *, t_mul, Th)	// 行列との積
				default:
					__assume(false)
			}
			#undef DEF_TEST
		}
	}
}
