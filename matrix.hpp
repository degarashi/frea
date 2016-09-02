#pragma once
#include "vector.hpp"
#include "angle.hpp"
#include "meta/compare.hpp"
#include "ieee754.hpp"

namespace frea {
	template <class VW, int M>
	using wrapM_t = wrapM_spec<VW, M, VW::size>;

	struct NoInverseMatrix : std::runtime_error {
		using std::runtime_error::runtime_error;
	};

	// ベクトル演算レジスタクラスをM方向に束ねたもの
	/*!
		\tparam	VW		各行のベクトル型
		\tparam	M		行数
		\tparam	S		本来の型
	*/
	template <class VW, int M, class S>
	class wrapM : public op::PlusMinus<S>, public op::Ne<S> {
		public:
			using op_t = op::PlusMinus<S>;
			constexpr static bool is_integral = VW::is_integral;
			constexpr static int dim_m = M,
								dim_n = VW::size,
								dim_min = Arithmetic<dim_m, dim_n>::less,
								bit_width = VW::bit_width;
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
			template <std::size_t... Idx, std::size_t... IdxE>
			static spec_t _Diagonal(const value_t& v, std::index_sequence<Idx...>, std::index_sequence<IdxE...>) {
				spec_t ret;
				_Dummy((ret.v[Idx].template initAt<Idx>(v), 0)...);
				_Dummy((ret.v[IdxE+dim_min] = vec_t::Zero())...);
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
				for(int i=0 ; i<dim_m ; i++) {
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
				using op_t::operator op;
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

			bool operator == (const wrapM& m) const {
				for(int i=0 ; i<dim_m ; i++) {
					if(v[i] != m.v[i])
						return false;
				}
				return true;
			}

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
				using ret_t = type_cn<ToM, ToN>;
				ret_t ret;
				for(int i=0 ; i<dim_m ; i++)
					ret.v[i] = v[i].template convert<ToN>();
				for(int i=dim_m ; i<ToM ; i++)
					ret.v[i] = ret_t::vec_t::Zero();
				return ret;
			}
			// 大きいサイズへの変換
			// 隙間をゼロで埋め、指定位置のみ任意の値を書き込む
			template <int ToM,
					 int ToN,
					 int Pos,
					 ENABLE_IF((ToM>dim_m))>
			auto convertI(const value_t& v) const {
				using ret_t = type_cn<ToM, ToN>;
				ret_t ret;
				for(int i=0 ; i<dim_m ; i++)
					ret.v[i] = this->v[i].template convertI<ToN, (Pos<dim_m ? Pos : 0)>(v);
				for(int i=dim_m ; i<ToM ; i++)
					ret.v[i] = ret_t::vec_t::Zero();
				// 指定位置の要素に値をセット
				if(Pos >= dim_m)
					ret.v[Pos].template setAt<Pos>(v);
				return ret;
			}
			static spec_t Identity() {
				return Diagonal(1);
			}
			static spec_t Diagonal(const value_t& v) {
				return _Diagonal(v,
					std::make_index_sequence<dim_min>(),
					std::make_index_sequence<dim_m-dim_min>()
				);
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
	struct wrapM_spec : wrapM<VW, M, wrapM_spec<VW,M,N>> {
		using base_t = wrapM<VW, M, wrapM_spec<VW,M,N>>;
		using base_t::base_t;

		wrapM_spec() = default;
	};
	#define DEF_FUNC(name, nameC) \
		void name() { \
			*this = nameC(); \
		}

	template <class R, int M, int N, bool A>
	using SMat_t = MatT_spec<SVec_t<R, N, A>, M, N>;
	template <class T, int M, int N>
	using RMat_t = MatT_spec<RVec_t<T, N>, M, N>;
	template <class T, int M, int N, bool A>
	using Mat_t = MatT_spec<Vec_t<T, N, A>, M, N>;

	//! 正方行列のみのメンバ関数を定義
	template <class VW, int S>
	class wrapM_spec<VW,S,S> : public wrapM<VW, S, wrapM_spec<VW,S,S>> {
		private:
			using base_t = wrapM<VW, S, wrapM_spec<VW,S,S>>;
			using this_t = wrapM_spec;
			using value_t = typename base_t::value_t;
			using mat_t = Mat_t<value_t, base_t::dim_m, base_t::dim_n, true>;
			// InfoにTranspose関数が用意されていればtrue
			template <class VW2,
					 std::size_t... Idx,
					 class B,
					 class = decltype(VW2::I::Transpose(std::declval<B>().v[Idx].m...))>
			static std::true_type hasmem(std::index_sequence<Idx...>);
			// そうでなければfalse
			template <class VW2, class B>
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
				const auto dummy = [](auto...) {};
				dummy((base_t::v[Idx].template store<false>(tmp[Idx], IConst<base_t::dim_n-1>()), 0)...);
				dummy(([&tmp](const int i){
					for(int j=i+1 ; j<base_t::dim_n ; j++) {
						std::swap(tmp[i][j], tmp[j][i]);
					}
				}(Idx), 0)...);
				return this_t((const typename base_t::value_t*)tmp, std::false_type());
			}
		public:
			using vec_t = typename base_t::vec_t;
			using base_t::base_t;
			wrapM_spec() = default;
	
			auto transposition() const& {
				const auto idx = std::make_index_sequence<S>();
				return _transposition(decltype(this->hasmem<VW, base_t>(idx))(), idx);
			}
			// ----------- アダプタ関数群 -----------
			value_t calcDeterminant() const { return mat_t(*this).calcDeterminant(); }
			this_t inversion() const { return mat_t(*this).inversion(); }
			this_t inversion(const value_t& det) const { return mat_t(*this).inversion(det); }
			DEF_FUNC(inverse, inversion)
			DEF_FUNC(transpose, transposition)
	};
	#undef DEF_FUNC

	template <bool C, class T>
	using ConstIf = std::conditional_t<C, std::add_const_t<T>, T>;

	template <class D, bool C>
	class ItrM :
		public std::iterator<
			std::input_iterator_tag,
			ConstIf<C, typename D::value_t>
		>
	{
		private:
			using base_t = std::iterator<std::input_iterator_tag, ConstIf<C, typename D::value_t>>;
			using D_t = ConstIf<C, D>;
			D_t&	_target;
			int		_cursor;
			constexpr static int dim_m = D::dim_m;
		public:
			ItrM(D_t& t, const int cur):
				_target(t),
				_cursor(cur)
			{}
			typename base_t::reference operator * () const {
				return _target.m[_cursor/dim_m][_cursor%dim_m];
			}
			ItrM& operator ++ () {
				++_cursor;
				return *this;
			}
			ItrM operator ++ (int) {
				const int cur = _cursor++;
				return ItrM(_target, cur);
			}
			typename base_t::pointer operator -> () const {
				return &(*this);
			}
			bool operator == (const ItrM& itr) const {
				return _cursor == itr._cursor;
			}
			bool operator != (const ItrM& itr) const {
				return _cursor != itr._cursor;
			}
	};

	template <class V, int M>
	class DataM {
		public:
			using value_t = typename V::value_t;
			using vec_t = V;
			constexpr static int dim_m = M,
								dim_n = vec_t::size,
								bit_width = vec_t::bit_width;
			constexpr static bool align = vec_t::align,
								is_integral = vec_t::is_integral;
			using column_t = typename vec_t::template type_cn<dim_m>;
		private:
			void _init(vec_t*) {}
			template <class... Ts>
			void _init(vec_t* dst, const vec_t& v, const Ts&... ts) {
				*dst = v;
				_init(dst+1, ts...);
			}
			template <std::size_t... Idx>
			void _initA(std::index_sequence<Idx...>, const value_t *src) {
				static_assert(sizeof...(Idx) == dim_n, "");
				for(int i=0 ; i<dim_m ; i++) {
					m[i] = vec_t(src[Idx]...);
					src += dim_n;
				}
			}
		public:
			vec_t	m[dim_m];
			DataM() = default;

			// --- iterator interface ---
			auto begin() { return ItrM<DataM, false>(*this, 0); }
			auto end() { return ItrM<DataM, false>(*this, dim_m*dim_n); }
			auto begin() const { return ItrM<DataM, true>(*this, 0); }
			auto end() const { return ItrM<DataM, true>(*this, dim_m*dim_n); }
			auto cbegin() const { return ItrM<DataM, true>(*this, 0); }
			auto cend() const { return ItrM<DataM, true>(*this, dim_m*dim_n); }

			vec_t& operator [](const int n) {
				return m[n];
			}
			const vec_t& operator [](const int n) const {
				return m[n];
			}
			// vec_tの配列[dim_m]で初期化
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
				vec_t tmp(t0);
				for(auto& v : m)
					v = tmp;
			}
			// from 内部形式
			DataM(const wrapM_spec<typename V::wrap_t, M, V::size>& w) {
				w.store(m);
			}
	};

	#define AsI(t) wrap_t(reinterpret_cast<const value_t*>((t).m), BConst<align>())
	template <class V, int M, class S>
	struct MatT : DataM<V,M>, op::Operator_Ne<S> {
		using op_t = op::Operator_Ne<S>;
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
		using type_cn = Mat_t<value_t, M2, N2, vec_t::template type_cn<N2>::align>;

		MatT() = default;
		constexpr MatT(const base_t& b): base_t(b) {}
		MatT(const wrap_t& w) {
			for(int i=0 ; i<dim_m ; i++)
				base_t::m[i] = w.v[i];
		}

		//! 指定した行と列を省いた物を出力
		auto cutRC(const int row, const int clm) const {
			type_cn<dim_m-1, dim_n-1>	ret;
			constexpr int width = dim_n,
						height = dim_m;
			// 左上
			for(int i=0 ; i<row ; i++) {
				for(int j=0 ; j<clm ; j++)
					ret.m[i][j] = this->m[i][j];
			}
			// 右上
			for(int i=0 ; i<row ; i++) {
				for(int j=clm+1 ; j<width ; j++)
					ret.m[i][j-1] = this->m[i][j];
			}
			// 左下
			for(int i=row+1 ; i<height ; i++) {
				for(int j=0 ; j<clm ; j++)
					ret.m[i-1][j] = this->m[i][j];
			}
			// 右下
			for(int i=row+1 ; i<height ; i++) {
				for(int j=clm+1 ; j<width ; j++)
					ret.m[i-1][j-1] = this->m[i][j];
			}
			return ret;
		}

		// ---------- < 行の基本操作 > ----------
		//! 行を入れ替える
		void rowSwap(const int r0, const int r1) {
			std::swap(this->m[r0], this->m[r1]);
		}
		//! ある行を定数倍する
		void rowMul(const int r0, const value_t& s) {
			this->m[r0] *= s;
		}
		//! ある行を定数倍した物を別の行へ足す
		void rowMulAdd(const int r0, const value_t& s, const int r1) {
			this->m[r1] += this->m[r0] * s;
		}
		// ---------- < 列の基本操作 > ----------
		//! 列を入れ替える
		void clmSwap(const int c0, const int c1) {
			for(int i=0 ; i<dim_m ; i++)
				std::swap(this->m[i][c0], this->m[i][c1]);
		}
		//! ある列を定数倍する
		void clmMul(const int c0, const value_t& s) {
			for(int i=0 ; i<dim_m ; i++)
				this->m[i][c0] *= s;
		}
		//! ある列を定数倍した物を別の行へ足す
		void clmMulAdd(const int c0, const value_t& s, const int c1) {
			for(int i=0 ; i<dim_m ; i++)
				this->m[i][c1] += this->m[i][c0] * s;
		}

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
		template <int At>
		column_t getColumn() const {
			static_assert(At < dim_n, "invalid position");
			column_t ret;
			for(int i=0 ; i<dim_m ; i++)
				ret[i] = this->m[i][At];
			return ret;
		}
		template <int At>
		void setColumn(const column_t& c) {
			static_assert(At < dim_n, "invalid position");
			for(int i=0 ; i<dim_m ; i++)
				this->m[i][At] = c[i];
		}

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
			using op_t::operator op;
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
	// 正方行列のみのメンバ関数を定義
	template <class V, int N, class S>
	class MatT_dspec<V,N,N,S> : public MatT<V,N,S> {
		public:
			using base_t = MatT<V,N,S>;
			using base_t::base_t;
			using spec_t = S;
			using value_t = typename V::value_t;
			using this_t = MatT_dspec;
		private:
			value_t _calcDeterminant(IConst<2>) const {
				// 公式で計算
				const Mat_t<value_t, 2, 2, false> m(*this);
				return m[0][0]*m[1][1] - m[0][1]*m[1][0];
			}
			template <int N2, ENABLE_IF((N2>2))>
			value_t _calcDeterminant(IConst<N2>) const {
				value_t res = 0,
						s = 1;
				// 部分行列を使って計算
				for(int i=0 ; i<N2 ; i++) {
					const auto mt = this->cutRC(0,i);
					res += this->m[0][i] * mt.calcDeterminant() * s;
					s *= -1;
				}
				return res;
			}
			spec_t  _inversion(const value_t& di, IConst<2>) const {
				spec_t ret;
				ret.m[0][0] = this->m[1][1] * di;
				ret.m[1][0] = -this->m[1][0] * di;
				ret.m[0][1] = -this->m[0][1] * di;
				ret.m[1][1] = this->m[0][0] * di;
				return ret;
			}
			template <int N2, ENABLE_IF((N2>2))>
			spec_t _inversion(const value_t& di, IConst<N2>) const {
				spec_t ret;
				const value_t c_val[2] = {1,-1};
				for(int i=0 ; i<base_t::dim_m ; i++) {
					for(int j=0 ; j<base_t::dim_n ; j++) {
						auto in_mat = this->cutRC(i,j);
						const value_t in_di = in_mat.calcDeterminant();
						ret.m[j][i] = c_val[(i+j)&1] * in_di * di;
					}
				}
				return ret;
			}
			constexpr static auto ZeroThreshold = Threshold<value_t>(0.6, 1);

		public:
			value_t calcDeterminant() const {
				return _calcDeterminant(IConst<base_t::dim_m>());
			}
			spec_t inversion() const {
				return inversion(calcDeterminant());
			}
			spec_t inversion(const value_t& det) const {
				constexpr auto Th = ZeroThreshold;
				if(std::abs(det) < Th)
					throw NoInverseMatrix("");
				return _inversion(1/det, IConst<base_t::dim_m>());
			}
			void invert() {
				*this = inversion();
			}

			// -------- アダプタ関数群(wrapM) --------
			spec_t transposition() const { return this->asInternal().transposition(); }
			void transpose() { *this = transposition(); }
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
