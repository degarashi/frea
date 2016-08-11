#pragma once
#include "vector.hpp"
#include "angle.hpp"
#include "meta/compare.hpp"

namespace frea {
	template <class VW, int M>
	using wrapM_t = wrapM_spec<VW, M, VW::size>;

	// ベクトル演算レジスタクラスをM方向に束ねたもの
	/*!
		\tparam	VW		各行のベクトル型
		\tparam	M		行数
		\tparam	S		本来の型
	*/
	template <class VW, int M, class S>
	class wrapM {
		public:
			constexpr static bool is_integral = VW::is_integral;
			constexpr static int dim_m = M,
								dim_n = VW::size;
			using spec_t = S;
			using vec_t = VW;
			using value_t = typename vec_t::value_t;
			using column_t = typename vec_t::template type_cn<dim_m>;
			template <int M2, int N2>
			using type_cn = wrapM_t<typename vec_t::template type_cn<N2>, M2>;
		private:
			template <class WM>
			static void _MultipleLine(vec_t&, const vec_t&, const WM&, IConst<dim_n>) {}
			template <class WM, int N>
			static void _MultipleLine(vec_t& dst, const vec_t& v, const WM& m, IConst<N>) {
				auto tv = v;
				tv.template makeEquality<N>();
				dst += tv * m.v[N];
				_MultipleLine(dst, v, m, IConst<N+1>());
			}
			template <int At, std::size_t... Idx>
			auto _getColumn(std::index_sequence<Idx...>) const {
				return column_t((v[Idx].template pickAt<At>())...);
			}
			template <int At, std::size_t... Idx>
			void _setColumn(const column_t& c, std::index_sequence<Idx...>) {
				const auto dummy = [](auto&&...){};
				dummy(((v[Idx].template setAt<At>() = c.template pickAt<At>()), 0)...);
			}
			template <class... Ts>
			static void _Dummy(Ts&&...) {}
			template <std::size_t... Idx>
			static spec_t _Diagonal(const value_t& v, std::index_sequence<Idx...>) {
				spec_t ret;
				_Dummy((ret.v[Idx].template initAt<Idx>(v), 0)...);
				return ret;
			}
			template <class VW2, int N2, class S2>
			auto _mul(const wrapM<VW2,N2,S2>& m, std::true_type) const {
				using WM = std::decay_t<decltype(m)>;
				static_assert(WM::dim_m == dim_n, "");
				wrapM_t<VW2, dim_m> ret;
				value_t other[dim_m][VW2::size];
				for(int i=0 ; i<VW2::size ; i++)
					m.v[i].template store<false>(other[i], IConst<VW2::size-1>());
				for(int i=0 ; i<WM::dim_m ; i++) {
					value_t result[WM::dim_n] = {},
							ths[dim_n];
					v[i].template store<false>(ths, IConst<dim_n-1>());
					for(int j=0 ; j<WM::dim_n ; j++) {
						auto& dst = result[j];
						for(int k=0 ; k<dim_n ; k++)
							dst += ths[k] * other[k][j];
					}
					ret.v[i] = VW2(result, std::false_type());
				}
				return ret;
			}
			template <class VW2, int N2, class S2>
			auto _mul(const wrapM<VW2,N2,S2>& m, std::false_type) const {
				using WM = std::decay_t<decltype(m)>;
				static_assert(WM::dim_m == dim_n, "");
				wrapM_t<VW2, dim_m> ret;
				for(int i=0 ; i<WM::dim_m ; i++) {
					_MultipleLine(ret.v[i] = vec_t::Zero(), v[i], m, IConst<0>());
				}
				return ret;
			}
			template <std::size_t... Idx, ENABLE_IF(sizeof...(Idx) == dim_m)>
			column_t _mul_vecR(std::index_sequence<Idx...>, const vec_t& vc) const {
				return column_t(v[Idx].dot(vc)...);
			}
			template <std::size_t... Idx, ENABLE_IF(sizeof...(Idx) == dim_n)>
			vec_t _mul_vecL(std::index_sequence<Idx...>, const column_t& vc) const {
				vec_t ret = vec_t::Zero();
				_MultipleLine(ret, vc, *this, IConst<0>());
				return ret;
			}
		public:
			vec_t	v[dim_m];

			wrapM() = default;
			template <bool A>
			wrapM(const value_t* src, BConst<A>) {
				for(int i=0 ; i<dim_n ; i++) {
					v[i] = vec_t(src, BConst<A>());
					src += A ? vec_t::capacity : vec_t::size;
				}
			}

			#define DEF_OP(op) \
				template <class T, \
							ENABLE_IF((HasMethod_asInternal_t<T>{}))> \
				spec_t operator op (const T& t) const { \
					return *this op t.asInternal(); \
				} \
				template <class VW2, class S2> \
				spec_t operator op (const wrapM<VW2,M,S2>& m) const { \
					spec_t ret; \
					for(int i=0 ; i<dim_m ; i++) \
						ret.v[i] = v[i] op m.v[i]; \
					return ret; \
				} \
				spec_t operator op (const value_t& t) const { \
					spec_t ret; \
					const vec_t tmp(t); \
					for(int i=0 ; i<dim_m ; i++) \
						ret.v[i] = v[i] op tmp; \
					return ret; \
				} \
				template <class T> \
				spec_t& operator op##= (const T& t) & { \
					return static_cast<spec_t&>(*this = *this op t); \
				}
			DEF_OP(+)
			DEF_OP(-)
			#undef DEF_OP
			template <class VW2, class S2>
			spec_t& operator = (const wrapM<VW2,dim_m,S2>& m) {
				for(int i=0 ; i<dim_m ; i++)
					v[i] = m.v[i];
				return *this;
			}
			auto operator * (const value_t& t) const {
				spec_t ret;
				vec_t tmp(t);
				for(int i=0 ; i<dim_m ; i++)
					ret.v[i] = v[i] * tmp;
				return ret;
			}
			template <class T = value_t,
					 ENABLE_IF(std::is_floating_point<T>{})>
			auto operator / (const value_t& t) const {
				return *this * (1 / t);
			}
			template <class T = value_t,
					 ENABLE_IF(std::is_integral<T>{})>
			auto operator / (const value_t& t) const {
				spec_t ret;
				const vec_t tmp(t);
				for(int i=0 ; i<dim_m ; i++)
					ret.v[i] = v[i] / tmp;
				return ret;
			}
			template <class VW2, int N2, class S2>
			auto operator * (const wrapM<VW2,N2,S2>& m) const {
				return _mul(m, IsTuple_t<VW2>());
			}
			template <class T>
			spec_t& operator *= (const T& m) {
				return static_cast<spec_t&>(*this = *this * m);
			}
			// 右からベクトルを掛ける
			column_t operator * (const vec_t& vc) const {
				return _mul_vecR(std::make_index_sequence<dim_m>(), vc);
			}
			// 左からベクトルを掛ける
			vec_t _pre_mul(const column_t& vc) const {
				return _mul_vecL(std::make_index_sequence<dim_n>(), vc);
			}

			bool operator == (const wrapM& m) const;
			bool operator != (const wrapM& m) const;

			template <class VD>
			void store(VD* dst) const {
				for(auto& t : v) {
					t.template store<VD::align>(dst->m, IConst<dim_n-1>());
					++dst;
				}
			}
			// 小さいサイズへの変換
			template <int ToM,
					 int ToN,
					 ENABLE_IF((ToM<=dim_m))>
			decltype(auto) convert() const {
				return reinterpret_cast<const type_cn<ToM, ToN>&>(*this);
			}
			template <int ToM,
					 int ToN,
					 int Pos,
					 ENABLE_IF((ToM<=dim_m))>
			decltype(auto) convertI(const value_t&) const {
				return reinterpret_cast<const type_cn<ToM, ToN>&>(*this);
			}
			// 大きいサイズへの変換
			// 隙間をゼロで埋める
			template <int ToM,
					 int ToN,
					 ENABLE_IF((ToM>dim_m))>
			auto convert() const {
				type_cn<ToM, ToN> ret;
				for(int i=0 ; i<dim_m ; i++)
					ret.v[i] = v[i].template convert<ToN>();
				for(int i=dim_m ; i<ToM ; i++)
					ret.v[i] = vec_t::I::Zero();
				return ret;
			}
			// 大きいサイズへの変換
			// 隙間をゼロで埋め、指定位置のみ任意の値を書き込む
			template <int ToM,
					 int ToN,
					 int Pos,
					 ENABLE_IF((ToM>dim_m))>
			auto convertI(const value_t& v) const {
				type_cn<ToM, ToN> ret;
				for(int i=0 ; i<dim_m ; i++)
					ret.v[i] = this->v[i].template convertI<ToN, (Pos<dim_m ? Pos : 0)>(v);
				for(int i=dim_m ; i<ToM ; i++)
					ret.v[i] = vec_t::I::Zero();
				// 指定位置の要素に値をセット
				if(Pos >= dim_m)
					ret.v[Pos].template setAt<Pos>(v);
				return ret;
			}
			static spec_t Identity() {
				return Diagonal(1);
			}
			static spec_t Diagonal(const value_t& v) {
				return _Diagonal(v, std::make_index_sequence<dim_m>());
			}
			template <int At>
			const auto& getRow() const {
				static_assert(At < dim_m, "invalid position");
				return v[At];
			}
			template <int At>
			auto& getRow() {
				static_assert(At < dim_m, "invalid position");
				return v[At];
			}
			template <int At>
			auto getColumn() const {
				static_assert(At < dim_n, "invalid position");
				return _getColumn<At>(std::make_index_sequence<dim_m>());
			}
			template <int At>
			void setColumn(const column_t& c) {
				static_assert(At < dim_n, "invalid position");
				_setColumn<At>(c, std::make_index_sequence<dim_m>());
			}
	};

	template <class VW, int M, int N>
	struct wrapM_spec : wrapM<VW, N, wrapM_spec<VW,M,N>> {
		using base_t = wrapM<VW, N, wrapM_spec<VW,M,N>>;
		using base_t::base_t;

		wrapM_spec() = default;
	};
	#define DEF_FUNC(name, nameC) \
		void name() { \
			*this = nameC(); \
		}

	//! 正方行列のみのメンバ関数を定義
	template <class VW, int S>
	class wrapM_spec<VW,S,S> : public wrapM<VW, S, wrapM_spec<VW,S,S>> {
		private:
			using base_t = wrapM<VW, S, wrapM_spec<VW,S,S>>;
			using this_t = wrapM_spec;
			// InfoにTranspose関数が用意されていればtrue
			template <class I2,
					 std::size_t... Idx,
					 class B,
					 class = decltype(I2::Transpose(std::declval<B>().v[Idx].m...))>
			static std::true_type hasmem(std::index_sequence<Idx...>);
			// そうでなければfalse
			template <class I2, class B>
			static std::false_type hasmem(...);
			//! 効率の良い転置アルゴリズムが使用可能な場合はそれを使用
			template <std::size_t... Idx>
			auto _transposition(std::true_type, std::index_sequence<Idx...>) const {
				auto ret = *this;
				VW::I::Transpose((ret.v[Idx].m)...);
				return ret;
			}
			//! 効率の良い転置アルゴリズムが定義されていない場合は地道に要素を交換する
			template <std::size_t... Idx>
			auto _transposition(std::false_type, std::index_sequence<Idx...>) const {
				typename base_t::value_t tmp[S][S];
				auto dummy = [](auto...) {};
				dummy((base_t::v[Idx].template store<false>(tmp[Idx], IConst<base_t::dim_n>()), 0)...);
				dummy(([&tmp](const int i){
					for(int j=0 ; j<(base_t::dim_n-i) ; j++) {
						std::swap(tmp[i][j], tmp[j][i]);
					}
				}(Idx), 0)...);
				return this_t((const typename base_t::value_t*)tmp, std::false_type());
			}
		public:
			using base_t::base_t;
			wrapM_spec() = default;
	
			auto transposition() const& {
				const auto idx = std::make_index_sequence<S>();
				return _transposition(decltype(this->hasmem<VW::I, base_t>(idx))(), idx);
			}
			DEF_FUNC(transpose, transposition)
	};
	#undef DEF_FUNC

	template <class V, int M>
	class DataM {
		public:
			using value_t = typename V::value_t;
			using vec_t = V;
			constexpr static int dim_m = M,
								dim_n = vec_t::size;
			constexpr static bool align = vec_t::align;
			using column_t = typename vec_t::template type_cn<dim_m>;
		private:
			using VData = typename V::base_t;
			void _init(VData*) {}
			template <class... Ts>
			void _init(VData* dst, const VData& v, const Ts&... ts) {
				*dst = v;
				_init(dst+1, ts...);
			}
			template <std::size_t... Idx>
			void _initA(std::index_sequence<Idx...>, const value_t *src) {
				static_assert(sizeof...(Idx) == dim_n, "");
				for(int i=0 ; i<dim_m ; i++) {
					m[i] = VData(src[Idx]...);
					src += dim_n;
				}
			}
		public:
			VData	m[dim_m];
			DataM() = default;

			auto& operator [](const int n) {
				return m[n];
			}
			const auto& operator [](const int n) const {
				return m[n];
			}
			// VDataの配列[dim_m]で初期化
			template <class... Ts,
					ENABLE_IF((sizeof...(Ts)==dim_m))>
			DataM(const Ts&... ts) {
				_init(m, ts...);
			}
			// value_tの配列で初期化
			template <class... Ts,
					ENABLE_IF((sizeof...(Ts)==dim_m*dim_n))>
			DataM(const Ts&... ts) {
				const value_t tmp[dim_m*dim_n] = {static_cast<value_t>(ts)...};
				_initA(std::make_index_sequence<dim_n>(), tmp);
			}
			// 全て同一のvalue_tで初期化
			explicit DataM(const value_t& t0) {
				VData tmp(t0);
				for(auto& v : m)
					v = tmp;
			}
			// from 内部形式
			DataM(const wrapM_spec<typename V::wrap_t, M, V::size>& w) {
				w.store(m);
			}
	};

	template <class R, int M, int N, bool A>
	using SMat_t = MatT_spec<SVec_t<R, N, A>, M, N>;
	template <class T, int M, int N>
	using RMat_t = MatT_spec<RVec_t<T, N>, M, N>;
	template <class T, int M, int N, bool A>
	using Mat_t = MatT_spec<Vec_t<T, N, A>, M, N>;

	#define AsI(t) wrap_t(reinterpret_cast<const value_t*>((t).m), BConst<align>())
	template <class V, int M, class S>
	struct MatT : DataM<V,M> {
		using spec_t = S;
		using base_t = DataM<V,M>;
		using base_t::base_t;
		using vec_t = V;
		using wrap_t = wrapM_t<typename V::wrap_t, M>;
		using value_t = typename V::value_t;
		constexpr static int dim_m = base_t::dim_m,
							dim_n = base_t::dim_n,
							dim_min = Arithmetic<dim_m,dim_n>::less;
		constexpr static int size = dim_m,
							lower_size = dim_n;
		constexpr static bool align = base_t::align;
		using column_t = typename V::template type_cn<dim_m>;
		using vec_min = typename vec_t::template type_cn<dim_min>;
		template <int M2, int N2>
		using type_cn = Mat_t<typename vec_t::template type_cn<N2>, M2, N2,
							vec_t::template type_cn<N2>::align>;

		MatT() = default;
		constexpr MatT(const base_t& b): base_t(b) {}
		MatT(const wrap_t& w) {
			for(int i=0 ; i<dim_m ; i++)
				base_t::m[i] = w.v[i];
		}

		//! 指定した行と列を省いた物を出力
		type_cn<dim_m-1, dim_n-1> cutRC(int row, int clm) const;

		//! 行列との積算 (右から掛ける)
		/*! 列ベクトルとして扱う = ベクトルを転置して左から行ベクトルを掛ける */
		vec_t operator * (const vec_t& v) const;
		// ---------- < 行の基本操作 > ----------
		//! 行を入れ替える
		void rowSwap(int r0, int r1);
		//! ある行を定数倍する
		void rowMul(int r0, float s);
		//! ある行を定数倍した物を別の行へ足す
		void rowMulAdd(const int r0, const value_t& s, const int r1);
		// ---------- < 列の基本操作 > ----------
		//! 列を入れ替える
		void clmSwap(int c0, int c1);
		//! ある列を定数倍する
		void clmMul(int c0, float s);
		//! ある列を定数倍した物を別の行へ足す
		void clmMulAdd(int c0, float s, int c1);

		//! ある行の要素が全てゼロか判定 (誤差=EPSILON込み)
		bool isZeroRowEps(int n) const;
		bool isZeroRow(int n) const;
		//! 零行列か判定
		bool isZero() const;
		//! 各行を正規化する (最大の係数が1になるように)
		void rowNormalize();
		//! 被約形かどうか判定
		bool isEchelon() const;
		//! 被約形にする
		/*! \return 0の行数 */
		int rowReduce();

		static MatT Diagonal(const value_t& v) {
			return wrap_t::Diagonal(v);
		}
		template <class... Ts, ENABLE_IF((sizeof...(Ts)==dim_min))>
		static MatT Diagonal(const value_t& v0, const Ts&... ts);
		static MatT Identity() {
			return wrap_t::Diagonal(1);
		}

		template <int At>
		const auto& getRow() const {
			static_assert(At < dim_m, "invalid position");
			return this->m[At];
		}
		template <int At>
		auto& getRow() {
			static_assert(At < dim_m, "invalid position");
			return this->m[At];
		}
		// template <int At>
		// column_t getColumn() const {
		// 	static_assert(At < dim_n, "invalid position");
		// 	return {this->m[Idx][At]...};
		// }
		// template <int At>
		// void setColumn(const column_t& c) {
		// 	static_assert(At < dim_n, "invalid position");
		// 	dummy(((this->m[Idx][At] = c.data[At]), 0)...);
		// }

		constexpr operator wrap_t() const {
			return AsI(*this);
		}
		wrap_t asInternal() const {
			return AsI(*this);
		}
		#define DEF_OP(op) \
			auto operator op (const wrap_t& m) const { \
				return AsI(*this) op m; \
			} \
			auto operator op (const value_t& v) const { \
				return AsI(*this) op v; \
			} \
			template <class T> \
			MatT& operator op##= (const T& m) { \
				return *this = *this op m; \
			}
		DEF_OP(+)
		DEF_OP(-)
		DEF_OP(*)
		DEF_OP(/)
		DEF_OP(&)
		DEF_OP(|)
		DEF_OP(^)
		#undef DEF_OP
		bool operator == (const MatT& m) const {
			return AsI(*this) == m;
		}
		bool operator != (const MatT& m) const {
			return !(this->operator == (m));
		}

		// ベクトルを左から掛ける
		auto _pre_mul(const column_t& vc) const {
			return vc * AsI(*this);
		}
		static spec_t Translation(const vec_min& v);
		static spec_t Scaling(const vec_min& v);

		//! サイズ変換
		//! 足りない要素はゼロで埋める
		template <int M2, int N2>
		auto convert() const {
			return AsI(*this).template convert<M2,N2>();
		}
		//! サイズ変換
		//! 足りない要素はゼロで埋め、対角線上要素を任意の値で初期化
		template <int M2, int N2, int Pos>
		auto convertI(const value_t& vi) const {
			return AsI(*this).template convertI<M2,N2,Pos>(vi);
		}
	};
	#undef AsI
	template <class V, int M, int N, class S>
	struct MatT_dspec : MatT<V,M,S> {
		using base_t = MatT<V,M,S>;
		using base_t::base_t;
	};
	template <class V, int N, class S>
	struct MatT_dspec<V,N,N,S> : MatT<V,N,S> {
		using base_t = MatT<V,N,S>;
		using base_t::base_t;
		using spec_t = S;
		using value_t = typename V::value_t;

		spec_t transposition() const { return this->asInternal().transposition(); }
		spec_t transpose() { return this->asInternal().transpose(); }
		value_t calcDeterminant() const;
		bool inversion(base_t& dst) const;
		bool inversion(base_t& dst, const value_t& det) const;
		bool invert();
	};
	template <class V, int M, int N>
	struct MatT_spec : MatT_dspec<V,M,N, MatT_spec<V,M,N>> {
		using base_t = MatT_dspec<V,M,N, MatT_spec<V,M,N>>;
		using base_t::base_t;
	};
}
#include "include/mat_d2.hpp"
#include "include/mat_d3.hpp"
#include "include/mat_d4.hpp"
namespace frea {
	#if SSE >= 2
		#define DEF_RM(n) \
			using RMat##n = RMat_t<float, n, n>;
		DEF_RM(2)
		DEF_RM(3)
		DEF_RM(4)
		#undef DEF_RM
	#endif

	#define DEF_M(n) \
		using Mat##n = Mat_t<float, n, n, false>; \
		using IMat##n = Mat_t<int32_t, n, n, false>; \
		using DMat##n = Mat_t<double, n, n, false>; \
		using AMat##n = Mat_t<float, n, n, true>;
	DEF_M(2)
	DEF_M(3)
	DEF_M(4)
	#undef DEF_M
}
