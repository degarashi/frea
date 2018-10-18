#pragma once
#include "lubee/src/error.hpp"
#include "lubee/src/meta/enable_if.hpp"
#include "detect_type.hpp"
#include "lubee/src/meta/index_sequence.hpp"
#include <type_traits>
#include <algorithm>
#include "lubee/src/meta/check_macro.hpp"
#include "lubee/src/meta/boolean.hpp"
#include "lubee/src/meta/compare.hpp"
#include "lubee/src/operators.hpp"
#include "lubee/src/tostring.hpp"

DEF_HASMETHOD(asInternal)

namespace frea {
	//! 演算レジスタに関する情報
	template <class R>
	struct info;

	namespace detail {
		// レジスタが要素数を内包できればそれを、そうでなければTupleに収める
		template <class R,
				 int N,
				 ENABLE_IF((info<R>::capacity < N))>
		auto DetectTup(const R&, lubee::IConst<N>) -> tup_spec<wrap_spec<R,info<R>::capacity>, N>;
		template <class R,
				 int N,
				 ENABLE_IF((info<R>::capacity >= N))>
		auto DetectTup(const R&, lubee::IConst<N>) -> wrap_spec<R,N>;
		template <class T, int N, bool A, int N2>
		auto DetectTup(const Data<T,N,A>&, lubee::IConst<N2>) -> wrap_spec<Data<T,N2,A>,N2>;
	}
	template <class R, int N>
	using Wrap_t = decltype(detail::DetectTup(std::declval<R>(), lubee::IConst<N>()));

	//! 演算レジスタラッパー
	/*!
		\tparam R	レジスタ型
		\tparam D	要素数
		\tparam S	本来の型
	*/
	template <class R, int D, class S>
	struct wrap : lubee::op::Operator_Ne<S> {
		using op_t = lubee::op::Operator_Ne<S>;
		using spec_t = S;
		using reg_t = R;
		using I = info<R>;
		using value_t = typename I::value_t;
		template <bool A>
		using data_t = Data<value_t, D, A>;
		constexpr static bool is_integral = static_cast<bool>(std::is_integral<value_t>{});
		//! 格納する予定の要素数
		constexpr static int size = D;
		//! 演算レジスタが格納できる要素数
		constexpr static int capacity = I::capacity,
							bit_width = sizeof(value_t)*8;
		//! 要素数の読み替え
		template <int D2>
		using type_cn = Wrap_t<reg_t, D2>;
		//! 内部表現値の読み替え
		template <class R2>
		using reg_cn = wrap_spec<R2, size>;

		//! コンパイル時チェック: 要素数がレジスタのサイズを超えていないか
		using Chk_Size = std::enable_if_t<(size<=capacity && size>0)>*;

		//! 実際に値を保持する型
		reg_t	 m;

		wrap() = default;
		decltype(auto) asInternal() const noexcept { return *this; }
		operator const reg_t& () const noexcept { return m; }
		wrap(const reg_t& r) noexcept: m(r) {}
		// 違うRegisterTypeの時は変換をかける
		template <class R2, class S2>
		wrap(const wrap<R2, size, S2>& w) noexcept:
			m(I::Convert(w.m))
		{}
		// 数値を指定されたら全要素に代入
		explicit wrap(const value_t& v) noexcept:
			m(I::Set1(v))
		{}
		// 要素数に満たない引数による初期化 (残りはゼロで埋める)
		template <class... Ts,
				 ENABLE_IF((capacity-2-sizeof...(Ts) > 0))>
		wrap(const value_t& v0, const value_t& v1, const Ts&... ts) noexcept:
			wrap(std::make_index_sequence<capacity-2-sizeof...(Ts)>(), v0, v1, ts...)
		{}
		template <class... Ts,
				 std::size_t... Idx>
		wrap(std::index_sequence<Idx...>, const Ts&... ts) noexcept:
			wrap(ts..., ZeroValue<Idx>()...)
		{}
		template <std::size_t Id>
		static value_t ZeroValue() noexcept { return 0; }
		// 要素数ピッタリの引数による初期化
		template <class... Ts,
				 ENABLE_IF((sizeof...(Ts)==capacity))>
		wrap(const Ts&... ts) noexcept:
			m(I::SetR(ts...))
		{}
		// アラインメモリからの読み込み
		wrap(const value_t* src, std::true_type) noexcept {
			alignas(16) value_t tmp[size];
			for(int i=0 ; i<size ; i++)
				tmp[i] = *src++;
			m = I::Load(tmp);
		}
		// 非アラインメモリからの読み込み
		wrap(const value_t* src, std::false_type) noexcept {
			value_t tmp[size];
			for(int i=0 ; i<size ; i++)
				tmp[i] = *src++;
			m = I::LoadU(tmp);
		}
		// 各種演算定義
		#define DEF_OP(op, func) \
			template <class T, ENABLE_IF(HasMethod_asInternal_t<T>{})> \
			auto operator op (const T& t) const& noexcept { \
				return *this op t.asInternal(); } \
			spec_t operator op (const value_t& t) const& noexcept { \
				return *this op wrap(t); } \
			template <class R2, class S2> \
			spec_t operator op (const wrap<R2, size, S2>& w) const& noexcept { \
				return I::func(m, w.m); } \
			using op_t::operator op;
		DEF_OP(+, Add)
		DEF_OP(-, Sub)
		DEF_OP(*, Mul)
		DEF_OP(/, Div)
		DEF_OP(&, And)
		DEF_OP(|, Or)
		DEF_OP(^, Xor)
		#undef DEF_OP
		// 左から行列にベクトルを掛ける
		template <class M,
				 ENABLE_IF(is_wrapM<M>{})>
		auto operator * (const M& m) const& noexcept {
			return type_cn<M::dim_n>(m._pre_mul(*this));
		}

		bool operator == (const wrap& w) const noexcept {
			return I::Equal(maskH<size-1>(), w.maskH<size-1>());
		}
		//! 上位要素をマスク(=0)する
		template <int N>
		spec_t maskH() const noexcept {
			static_assert(N<size, "");
			return I::And(m, I::MaskH(lubee::IConst<N>()));
		}
		template <int N>
		static auto MaskL(lubee::IConst<N>) noexcept {
			return I::Xor(I::MaskH(lubee::IConst<N>()), I::One());
		}
		//! 下位要素をマスク(=0)する
		template <int N>
		spec_t maskL() const noexcept {
			static_assert(N<size, "");
			return I::And(m, MaskL(lubee::IConst<N>()));
		}

		//! 指定要素に任意の値をセット(他はいじらない)
		template <int N>
		void setAt(const value_t& v) noexcept {
			static_assert(N<size, "");
			const auto mask = I::PickAt(lubee::IConst<N>());
			const auto maskInv = I::Xor(I::One(), mask);
			m = I::Or(I::And(m, maskInv), I::And(mask, I::Set1(v)));
		}

		//! 指定要素へ引数の値をセット
		template <int N>
		void initAt(const value_t& v) noexcept {
			m = I::And(I::PickAt(lubee::IConst<N>()), I::Set1(v));
		}
		//! 指定要素の取得
		template <int N, ENABLE_IF((N<size))>
		auto pickAt() const noexcept {
			return I::template Pick<N>(m);
		}
		//! 指定要素で他の要素を埋める
		template <int N, ENABLE_IF((N<size))>
		void makeEquality() noexcept {
			m = I::Set1(I::template Pick<N>(m));
		}
		value_t dot(const wrap& w) const noexcept {
			return (*this*w).sumUp();
		}
		value_t average() const noexcept {
			return sumUp() / size;
		}
		value_t distance(const wrap& w) const noexcept {
			return std::sqrt(dist_sq(w));
		}
		//! 距離の二乗
		value_t dist_sq(const wrap& w) const noexcept {
			const auto tv = w - *this;
			return tv.len_sq();
		}
		//! 要素ごとの最小値
		spec_t getMin(const wrap& w) const noexcept {
			return I::Min(m, w.m);
		}
		//! 要素ごとに最小値で上書き
		void selectMin(const wrap& w) noexcept {
			*this = getMin(w);
		}
		//! 要素ごとの最大値
		spec_t getMax(const wrap& w) const noexcept {
			return I::Max(m, w.m);
		}
		//! 要素ごとに最大値で上書き
		void selectMax(const wrap& w) noexcept {
			*this = getMax(w);
		}

		spec_t operator - () const noexcept {
			return *this * wrap(-1);
		}
		/*! \return 正規化前のベクトル長 */
		value_t normalize() noexcept {
			const value_t len = length();
			*this /= len;
			return len;
		}
		spec_t normalization() const noexcept {
			return *this * I::Reciprocal(length());
		}
		value_t length() const noexcept {
			return std::sqrt(len_sq());
		}
		//! ベクトル長の二乗
		value_t len_sq() const noexcept {
			return dot(*this);
		}
		//! 要素に1つでもNaNが含まれているか
		bool isNaN() const noexcept {
			return I::IsNaN(maskH<size-1>());
		}
		//! 要素に1つでもNaNまたはinfが含まれているか
		bool isOutstanding() const noexcept {
			return I::IsOutstanding(maskH<size-1>());
		}
		spec_t saturation(const value_t& vMin, const value_t& vMax) const noexcept {
			auto tm = I::Max(m, I::Set1(vMin));
			tm = I::Min(tm, I::Set1(vMax));
			return tm;
		}
		//! 線形補間
		spec_t l_intp(const wrap& w, const value_t& r) const noexcept {
			return *this + (w - *this) * r;
		}

		//! 各要素を足し合わせる
		auto sumUp() const noexcept {
			return I::SumUp(maskH<size-1>());
		}
		//! メモリへの書き出し
		/*!
			\tparam A	アラインメント
			\tparam N	書き出す要素数(0=最下位要素のみ)
		*/
		template <bool A, int N>
		void store(value_t* dst, lubee::IConst<N>) const noexcept {
			static_assert(N<size, "");
			I::Store(dst, m, lubee::BConst<A>(), lubee::IConst<N>());
		}
		//! レジスタ型の読み替え
		template <class R2>
		auto cast() const noexcept {
			return reg_cn<R2>(I::Cast(m));
		}
		//! 小さいサイズへの変換
		template <int ToN,
				 ENABLE_IF((ToN<=size))>
		decltype(auto) convert() const noexcept {
			return reinterpret_cast<const type_cn<ToN>&>(*this);
		}
		template <int ToN,
				 int Pos,
				 ENABLE_IF((ToN<=size))>
		decltype(auto) convertI(const value_t&) const noexcept {
			return reinterpret_cast<const type_cn<ToN>&>(*this);
		}
		//! 大きいサイズへの変換
		/*! 隙間をゼロで埋める */
		template <int ToN,
				 ENABLE_IF((ToN>size))>
		auto convert() const noexcept {
			return type_cn<ToN>(maskH<size-1>());
		}
		//! 大きいサイズへの変換
		/*! 隙間をゼロで埋め、指定位置のみ任意の値を書き込む */
		template <int ToN,
				 int Pos,
				 ENABLE_IF((ToN>size))>
		auto convertI(const value_t& v) const noexcept {
			static_assert(Pos<ToN, "");
			auto ret = convert<ToN>();
			// 指定位置の要素に値をセット
			if(Pos >= size)
				ret.template setAt<Pos>(v);
			return ret;
		}
		//! 要素ごとの絶対値
		spec_t absolute() const noexcept {
			return I::Absolute(m);
		}
		//! 要素が全てゼロであるかの判定
		bool isZero(const value_t& th) const noexcept {
			auto r = I::Lt(absolute().template maskH<size-1>(), wrap(th)),
				full = I::One(),
				one = I::Set1(1);
			return wrap(I::And(one, I::Xor(r, full))).sumUp() == 0;
		}
		//! 要素の中で最小の値
		auto getMinValue() const noexcept {
			return I::GetMinValue(
				maskH<size-1>() |
				spec_t(I::And(I::Set1(std::numeric_limits<value_t>::max()), MaskL(lubee::IConst<size-1>())))
			);
		}
		//! 要素の中で最大の値
		auto getMaxValue() const noexcept {
			return I::GetMaxValue(
				maskH<size-1>() |
				spec_t(I::And(I::Set1(std::numeric_limits<value_t>::lowest()), MaskL(lubee::IConst<size-1>())))
			);
		}
		//! 要素の最大値で他を除算
		void linearNormalize() noexcept {
			*this = linearNormalization();
		}
		spec_t linearNormalization() const noexcept {
			const value_t mv = absolute().getMaxValue();
			return *this / wrap(mv);
		}
		static spec_t Zero() noexcept { return I::Zero(); }
	};
	// wrapクラスに要素数固有の関数などを付加
	template <class R, int N>
	struct wrap_spec : wrap<R,N, wrap_spec<R,N>> {
		using base_t = wrap<R,N, wrap_spec<R,N>>;
		using base_t::base_t;
		wrap_spec() = default;
		wrap_spec(const base_t& b): base_t(b) {}
	};
	// tupクラスに要素数固有の関数などを付加
	template <class T, int N>
	struct tup_spec : tup<T,N, tup_spec<T,N>> {
		using base_t = tup<T,N, tup_spec<T,N>>;
		using base_t::base_t;
		tup_spec() = default;
		tup_spec(const base_t& b): base_t(b) {}
	};
}
#include "include/wrap_d2.hpp"
#include "include/wrap_d3.hpp"
#include "include/wrap_d4.hpp"
namespace frea{
	//! 要素が1レジスタでは収まりきらない場合に配列化
	/*!
		\tparam W wrapクラス
		\tparam N 要素数
		\tparam S 元のクラス
	*/
	template <class W, int N, class S>
	struct tup : lubee::op::Operator_Ne<S> {
		using op_t = lubee::op::Operator_Ne<S>;
		using wrap_t = W;
		using spec_t = S;
		using I = typename wrap_t::I;
		using value_t = typename wrap_t::value_t;
		using reg_t = typename wrap_t::reg_t;
		constexpr static int size = N,										//!< 要素数
							w_capacity = wrap_t::capacity,					//!< レジスタ1つが保持できる要素数
							a_size = (size-1)/w_capacity+1,					//!< レジスタをいくつ保持しているか
							capacity = w_capacity * a_size,					//!< このクラスが保持できる最大要素数
							bit_width = sizeof(value_t)*8;
		constexpr static bool is_integral = wrap_t::is_integral;
		constexpr static int Mod = N % w_capacity,
							Rem = (Mod==0) ? w_capacity : Mod;
		using wrapTail_t = typename wrap_t::template type_cn<Rem>;

		auto getMaskedTail() const noexcept {
			return data[a_size-1].template maskH<Rem-1>();
		}
		auto& getTail() noexcept {
			return data[a_size-1];
		}
		const auto& getTail() const noexcept {
			return data[a_size-1];
		}

		wrap_t		data[a_size];

		//! 要素数ピッタリの値を指定し初期化
		template <
			class... Ts,
			ENABLE_IF((
				 lubee::meta::And<std::is_convertible<Ts,value_t>...>{} &&
					sizeof...(Ts)==capacity
			))
		>
		tup(const Ts&... ts) noexcept {
			// 一度配列に展開
			alignas(16) const value_t tmp[capacity] = {static_cast<value_t>(ts)...};
			load(tmp, std::true_type());
		}
		template <std::size_t>
		constexpr static value_t ZeroValue() noexcept { return 0; }
		template <std::size_t... Idx, class... Ts>
		tup(std::index_sequence<Idx...>, const Ts&... ts) noexcept:
			tup(ts..., ZeroValue<Idx>()...) {}

		//! 要素数に満たない値を指定しての初期化 (=残りはゼロ)
		template <class... Ts,
				 ENABLE_IF((capacity-2-sizeof...(Ts)>0))>
		tup(const value_t& t0, const value_t& t1, const Ts&... ts) noexcept:
			tup(std::make_index_sequence<capacity-2-sizeof...(Ts)>(), t0, t1, ts...) {}
		//! メモリ配列から値を読み込む
		template <bool A>
		tup(const value_t* src, lubee::BConst<A>) noexcept {
			load(src, lubee::BConst<A>());
		}

		//! 要素数の読み替え
		template <int D2>
		using type_cn = Wrap_t<reg_t, D2>;

		auto dot(const tup& t) const noexcept {
			auto sum = wrap_t::Zero();
			for(int i=0 ; i<a_size-1 ; i++)
				sum += data[i] * t.data[i];
			sum += getMaskedTail() * t.getMaskedTail();
			return sum.sumUp();
		}
		auto average() const noexcept {
			return sumUp() / size;
		}
		auto sumUp() const noexcept {
			auto sum = wrap_t::Zero();
			for(int i=0 ; i<a_size-1 ; i++)
				sum += data[i];
			sum += getMaskedTail();
			return sum.sumUp();
		}
		auto distance(const tup& t) const noexcept {
			return std::sqrt(dist_sq(t));
		}
		auto dist_sq(const tup& t) const noexcept {
			const auto tmp = t - *this;
			return tmp.len_sq();
		}
		auto len_sq() const noexcept {
			return dot(*this);
		}
		spec_t getMin(const tup& t) const noexcept {
			spec_t ret;
			for(int i=0 ; i<a_size ; i++)
				ret.data[i] = data[i].getMin(t.data[i]);
			return ret;
		}
		void selectMin(const tup& t) noexcept {
			for(int i=0 ; i<a_size ; i++)
				data[i].selectMin(t.data[i]);
		}
		spec_t getMax(const tup& t) const noexcept {
			spec_t ret;
			for(int i=0 ; i<a_size ; i++)
				ret.data[i] = data[i].getMax(t.data[i]);
			return ret;
		}
		void selectMax(const tup& t) noexcept {
			for(int i=0 ; i<a_size ; i++)
				data[i].selectMax(t.data[i]);
		}
		auto length() const noexcept {
			return std::sqrt(len_sq());
		}
		auto normalization() const noexcept {
			return *this / length();
		}
		auto normalize() noexcept {
			const auto len = length();
			*this /= len;
			return len;
		}
		spec_t saturation(const value_t& vMin, const value_t& vMax) const noexcept {
			spec_t ret;
			for(int i=0 ; i<a_size ; i++)
				ret.data[i] = data[i].saturation(vMin, vMax);
			return ret;
		}
		spec_t l_intp(const tup& w, const value_t& r) const noexcept {
			spec_t ret;
			for(int i=0 ; i<a_size ; i++)
				ret.data[i] = data[i].l_intp(w.data[i], r);
			return ret;
		}
		bool isNaN() const noexcept {
			for(int i=0 ; i<a_size-1 ; i++) {
				if(data[i].isNaN())
					return true;
			}
			return getMaskedTail().isNaN();
		}
		bool isOutstanding() const noexcept {
			for(int i=0 ; i<a_size-1 ; i++) {
				if(data[i].isOutstanding())
					return true;
			}
			return getMaskedTail().isOutstanding();
		}
		spec_t operator - () const noexcept {
			return *this * spec_t(-1);
		}

		tup() = default;
		explicit tup(const value_t& t) noexcept {
			for(auto& d : this->data)
				d = wrap_t(t);
		}
		//! メモリからの読み込み
		template <bool A>
		void load(const value_t* src, lubee::BConst<A>) noexcept {
			for(auto& d : this->data) {
				d = wrap_t(src, lubee::BConst<A>());
				src += w_capacity;
			}
		}
		//! 他wrap or tuple形式からの変換
		template <class W2, ENABLE_IF(is_wrap<W2>{} || IsTuple_t<W2>{})>
		explicit tup(const W2& w) noexcept {
			// 一度メモリへ展開
			alignas(16) value_t tmp[W2::capacity];
			w.template store<true>(tmp, lubee::IConst<W2::size-1>());
			const auto* src = tmp;
			for(auto& d : this->data) {
				d = wrap_t(src, std::true_type());
				src += w_capacity;
			}
		}
		//! メモリへの書き出し
		template <bool A>
		void store(value_t* dst, lubee::IConst<N-1>) const noexcept {
			for(int i=0 ; i<a_size-1 ; i++) {
				this->data[i].template store<A>(dst, lubee::IConst<w_capacity-1>());
				dst += w_capacity;
			}
			getTail().template store<A>(dst, lubee::IConst<Rem-1>());
		}
		bool operator == (const tup& t) const noexcept {
			for(int i=0 ; i<a_size-1 ; i++) {
				if(this->data[i] != t.data[i])
					return false;
			}
			return getMaskedTail() == t.getMaskedTail();
		}
		#define DEF_OP(op) \
			template <class T, ENABLE_IF(HasMethod_asInternal_t<T>{})> \
			auto operator op (const T& t) const& noexcept { \
				return *this op t.asInternal(); } \
			spec_t operator op (const value_t& t) const& noexcept { \
				spec_t ret; \
				for(int i=0 ; i<a_size ; i++) \
					ret.data[i] = this->data[i] op t; \
				return ret; \
			} \
			spec_t operator op (const tup& t) const& noexcept { \
				spec_t ret; \
				for(int i=0 ; i<a_size ; i++) \
					ret.data[i] = this->data[i] op t.data[i]; \
				return ret; \
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
		//! 左から行列にベクトルを掛ける
		template <class M,
				 ENABLE_IF(is_wrapM<M>{})>
		type_cn<M::dim_n> operator * (const M& m) const& noexcept {
			return m._pre_mul(static_cast<const spec_t&>(*this));
		}

		//! 小さいサイズへの変換
		template <int ToN,
				 ENABLE_IF((ToN<=size))>
		decltype(auto) convert() const noexcept {
			return reinterpret_cast<const tup_spec<wrap_t,ToN>&>(*this);
		}
		template <int ToN,
				 int Pos,
				 ENABLE_IF((ToN<=size))>
		decltype(auto) convertI(const value_t&) const noexcept {
			return reinterpret_cast<const tup_spec<wrap_t,ToN>&>(*this);
		}
		//! 大きいサイズへの変換
		template <int ToN,
				 ENABLE_IF((ToN>size))>
		auto convert() const noexcept {
			tup_spec<wrap_t, ToN> ret;
			for(int i=0 ; i<a_size ; i++)
				ret.data[i] = this->data[i];
			ret.data[a_size-1] = getMaskedTail();
			const auto zero = wrap_t::I::Zero();
			constexpr auto AS = decltype(ret)::a_size;
			for(int i=a_size ; i<AS ; i++)
				ret.data[i] = zero;
			return ret;
		}
		template <int ToN,
				 int Pos,
				 ENABLE_IF((ToN>size))>
		auto convertI(const value_t& vi) const noexcept {
			static_assert(Pos<ToN, "");
			auto ret = convert<ToN>();
			if(Pos >= size)
				ret.data[Pos/w_capacity].template setAt<Pos % w_capacity>(vi);
			return ret;
		}
		template <int N2>
		void makeEquality() noexcept {
			const auto val = pickAt<N2>();
			const wrap_t tval(val);
			for(auto& v: this->data)
				v = tval;
		}
		template <int Pos>
		auto pickAt() const noexcept {
			static_assert(Pos<size, "");
			return this->data[Pos/a_size].template pickAt<Pos%a_size>();
		}
		template <int Pos>
		void initAt(const value_t& val) noexcept {
			const auto zero = wrap_t::Zero();
			for(auto& a : this->data)
				a = zero;
			this->data[Pos/a_size].template initAt<Pos%w_capacity>(val);
		}
		template <int Pos>
		spec_t maskH() const noexcept {
			spec_t ret(*this);
			constexpr auto From = Pos / w_capacity,
							Mod = Pos % w_capacity;
			const auto zero = wrap_t::Zero();
			for(int i=From+1 ; i<a_size ; i++)
				ret.data[i] = zero;
			ret.data[From] = this->data[From].template maskH<Mod>();
			return ret;
		}
		template <int Pos>
		spec_t maskL() const noexcept {
			spec_t ret(*this);
			constexpr auto From = Pos / w_capacity,
							Mod = Pos % w_capacity;
			const auto zero = wrap_t::Zero();
			for(int i=From-1 ; i>=0 ; i--)
				ret.data[i] = zero;
			ret.data[From] = this->data[From].template maskL<Mod>();
			return ret;
		}
		template <int Pos>
		void setAt(const value_t& val) noexcept {
			constexpr auto At = Pos / w_capacity,
							Mod = Pos % w_capacity;
			this->data[At].template setAt<Mod>(val);
		}
		spec_t absolute() const noexcept {
			spec_t ret(*this);
			for(auto& t : ret.data)
				t = I::Absolute(t);
			return ret;
		}
		bool isZero(const value_t& th) const noexcept {
			const auto tmp = absolute();
			for(int i=0 ; i<a_size-1 ; i++) {
				if(!tmp.data[i].isZero(th))
					return false;
			}
			return tmp.getMaskedTail().isZero(th);
		}
		value_t getMinValue() const noexcept {
			value_t m = data[0].getMinValue();
			for(int i=1 ; i<a_size-1 ; i++)
				m = std::min(m, data[i].getMinValue());
			m = std::min(m, wrapTail_t(getTail()).getMinValue());
			return m;
		}
		value_t getMaxValue() const noexcept {
			value_t m = data[0].getMaxValue();
			for(int i=1 ; i<a_size-1 ; i++)
				m = std::max(m, data[i].getMaxValue());
			m = std::max(m, wrapTail_t(getTail()).getMaxValue());
			return m;
		}
		void linearNormalize() noexcept {
			*this = linearNormalization();
		}
		spec_t linearNormalization() const noexcept {
			const value_t mv = absolute().getMaxValue();
			return *this / spec_t(mv);
		}
		static spec_t Zero() noexcept { return spec_t(0); }
	};
	// Dataクラスに要素数固有の関数などを付加
	template <class T, int N>
	struct Data_spec {
		T	m[N];
	};
	#define DEF_DATA(N, ...) \
		template <class T> \
		struct Data_spec<T,N> { \
			union { \
				T	m[N]; \
				struct { \
					T __VA_ARGS__; \
				}; \
			}; \
			Data_spec() = default; \
			template <class... Ts, ENABLE_IF((lubee::meta::And<std::is_convertible<Ts,T>...>{} && sizeof...(Ts)>=1))> \
			constexpr Data_spec(const Ts&... ts) noexcept: m{static_cast<T>(ts)...} {} \
		};
	DEF_DATA(1, x)
	DEF_DATA(2, x,y)
	DEF_DATA(3, x,y,z)
	DEF_DATA(4, x,y,z,w)
	#undef DEF_DATA

	static struct Set_t{} TagSet;
	static struct Mask_t{} TagMask;

	//! メモリ上でのデータ表現
	/*!
		\tparam T		内部表現型
		\tparam N		要素数
		\tparam A		アラインメント(true=16byte)
	*/
	template <class T, int N, bool A>
	struct alignas(A ? 16 : alignof(T)) Data : Data_spec<T,N> {
		constexpr static int size = N,
							capacity = size,
							bit_width = sizeof(T)*8;
		constexpr static bool align = A,
							is_integral = static_cast<bool>(std::is_integral<T>{});
		using base_t = Data_spec<T,N>;
		using value_t = T;
		using reg_t = Data;
		using I = info<Data>;
		//! 要素数, アラインメントの読み替え
		template <int N2, bool A2=align>
		using type_cn = Data<value_t,N2,A2>;

		// --- iterator interface ---
		value_t* begin() noexcept { return this->m; }
		value_t* end() noexcept { return this->m+N; }
		const value_t* begin() const noexcept { return this->m; }
		const value_t* end() const noexcept { return this->m+N; }
		const value_t* cbegin() const noexcept { return this->m; }
		const value_t* cend() const noexcept { return this->m+N; }

		Data() = default;
		template <bool A2>
		constexpr Data(const value_t* src, lubee::BConst<A2>) noexcept {
			for(auto& dst : base_t::m)
				dst = *src++;
		}
		const value_t& operator [](const int n) const noexcept { return this->m[n]; }
		value_t& operator [](const int n) noexcept { return this->m[n]; }
		//! 複数要素での初期化
		template <class... Ts,
				 ENABLE_IF((lubee::meta::And<std::is_convertible<Ts,value_t>...>{} && sizeof...(Ts)>1))>
		constexpr Data(const Ts&... ts) noexcept: base_t{static_cast<value_t>(ts)...} {}
		//! 1要素での初期化(内部呼び出し用)
		template <std::size_t... Idx, class T2>
		constexpr Data(std::index_sequence<Idx...>, const T2* t0) noexcept:
			base_t{static_cast<value_t>(t0[Idx])...} {}
		//! 1要素での初期化
		explicit constexpr Data(const value_t& t0) noexcept:
			Data(lubee::seq::Repeat_t<0,size>(), &t0) {}
		//! 指定要素以前をt0, 以降はゼロで初期化
		template <int Pos, class T2>
		Data(Mask_t, lubee::IConst<Pos>, const T2& t0) noexcept {
			for(int i=0 ; i<=Pos ; i++)
				this->m[i] = t0;
			for(int i=Pos+1 ; i<size ; i++)
				this->m[i] = 0;
		}
		//! 指定要素のみt0, 他はゼロで初期化
		template <int Pos, class T2>
		Data(Set_t, lubee::IConst<Pos>, const T2& t0) noexcept: base_t{} {
			this->m[Pos] = t0;
		}
		//! 内部形式からの初期化(wrap or tuple)
		template <class W,
				 class=decltype(std::declval<W>().template store<A>(base_t::m, lubee::IConst<size-1>()))>
		Data(const W& w) noexcept {
			constexpr int S = lubee::Arithmetic<size, W::size>::less;
			w.template store<A>(base_t::m, lubee::IConst<S-1>());
			for(int i=S ; i<size ; i++)
				base_t::m[i] = 0;
		}
		template <bool A2, int N2>
		void store(value_t* dst, lubee::IConst<N2>) const noexcept {
			const auto* src = base_t::m;
			int n = N2+1;
			while(--n >= 0) {
				*dst++ = *src++;
			}
		}

		template <int S,
				 ENABLE_IF(S<8)>
		static int32_t ToInt(lubee::IConst<S>);
		template <int S,
				 ENABLE_IF(S==8)>
		static int64_t ToInt(lubee::IConst<S>);
		template <class Typ>
		using ToInt_t = decltype(ToInt(lubee::IConst<sizeof(Typ)>()));

		#define DEF_OP(op) \
			template <class T2, bool A2> \
			Data operator op (const Data<T2,N,A2>& d) const noexcept { \
				Data ret; \
				for(int i=0 ; i<size ; i++) \
					ret.m[i] = base_t::m[i] op d.m[i]; \
				return ret; \
			}
		DEF_OP(+)
		DEF_OP(-)
		DEF_OP(*)
		DEF_OP(/)
		#undef DEF_OP

		#define DEF_OP(op) \
			template <class T2, bool A2> \
			Data operator op (const Data<T2,N,A2>& d) const noexcept { \
				Data ret; \
				for(int i=0 ; i<size ; i++) \
					ret.m[i] = ToInt_t<T>(base_t::m[i]) op ToInt_t<T2>(d.m[i]); \
				return ret; \
			}
		DEF_OP(&)
		DEF_OP(|)
		DEF_OP(^)
		#undef DEF_OP

		template <class Ar, class T2, int N2, bool A2>
		friend void serialize(Ar&, Data<T2,N2,A2>&);

		std::ostream& print(std::ostream& os) const {
			os << '[';
			bool bF = true;
			for(int i=0 ; i<size ; i++) {
				if(!bF)
					os << ", ";
				os << this->m[i];
				bF = false;
			}
			return os << ']';
		}
	};
	template <class R, int N, bool A>
	using SVec_t = VecT_spec<Wrap_t<R,N>, Data<typename info<R>::value_t, N, A>, N>;
	template <class T, int N>
	using RVec_t = VecT_spec<wrap_spec<Data<T,N,false>,N>, Data<T,N,false>, N>;

	template <int N, bool A, class T>
	RVec_t<T,N> info_detect(T);
}
#if SSE >= 2
	#include "include/sse.hpp"
#endif
#include "include/raw.hpp"
namespace frea {
	template <class W, class D, class S>
	struct VecT : D, lubee::op::Operator_Ne<S> {
		using op_t = lubee::op::Operator_Ne<S>;
		using spec_t = S;
		using base_t = D;
		using base_t::base_t;
		using wrap_t = W;
		using reg_t = typename wrap_t::reg_t;
		using value_t = typename wrap_t::value_t;

		constexpr static int size = W::size;
		constexpr static bool align = D::align;
		//! 要素数の読み替え
		template <int N2, bool A2=align>
		using type_cn = VecT_spec<Wrap_t<reg_t, N2>, typename base_t::template type_cn<N2,A2>, N2>;
		using Chk_Size = std::enable_if_t<(size>0)>*;

		auto _asInternal(std::false_type) const noexcept {
			return wrap_t(base_t::m, lubee::BConst<align>());
		}
		decltype(auto) _asInternal(std::true_type) const noexcept {
			return static_cast<const D&>(*this);
		}
		decltype(auto) asInternal() const noexcept {
			return _asInternal(std::is_same<base_t, wrap_t>());
		}
		constexpr operator wrap_t() const noexcept {
			return _asInternal(std::is_same<base_t, wrap_t>());
		}
		VecT() = default;
		template <class TA, ENABLE_IF(size==1 && (std::is_same<TA,value_t>{}))>
		VecT(const TA& t): base_t(t) {}
		constexpr VecT(const base_t& b) noexcept: base_t(b) {}
		template <class V,
				 ENABLE_IF(is_vector<V>{} && V::size==size)>
		VecT(const V& v) noexcept {
			for(int i=0 ; i<size ; i++)
				this->m[i] = v.m[i];
		}
		// 内部で別のレジスタ形式を内包しているベクトルからの変換
		template <class T,
				 ENABLE_IF((is_wrap<T>{} || IsTuple_t<T>{}) && T::size==size)>
		VecT(const T& t) noexcept {
			alignas(16) typename T::value_t tmp[size];
			t.template store<true>(tmp, lubee::IConst<size-1>());
			for(int i=0 ; i<size ; i++)
				this->m[i] = tmp[i];
		}
		bool operator == (const wrap_t& v) const noexcept {
			return asInternal() == v;
		}
		// ------------------- アダプタ関数 -------------------
		// 内部形式との演算子
		#define DEF_OP(op) \
			template <class W2> \
			auto operator op (const W2& w) const& noexcept { \
				return asInternal() op w; \
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
		auto dot(const wrap_t& w) const noexcept { return asInternal().dot(w); }
		value_t average() const noexcept { return asInternal().average(); }
		value_t distance(const wrap_t& p) const noexcept { return asInternal().distance(p); }
		value_t dist_sq(const wrap_t& p) const noexcept { return asInternal().dist_sq(p); }
		auto getMin(const wrap_t& p) const noexcept { return asInternal().getMin(p); }
		void selectMin(const wrap_t& p) noexcept {
			auto tmp = asInternal();
			tmp.selectMin(p);
			*this = tmp;
		}
		auto getMax(const wrap_t& p) const noexcept { return asInternal().getMax(p); }
		void selectMax(const wrap_t& p) noexcept {
			auto tmp = asInternal();
			tmp.selectMax(p);
			*this = tmp;
		}

		spec_t operator - () const noexcept { return -asInternal(); }
		value_t normalize() noexcept {
			auto tmp = asInternal();
			const value_t ret = tmp.normalize();
			*this = tmp;
			return ret;
		}
		spec_t normalization() const noexcept { return asInternal().normalization(); }
		value_t length() const noexcept { return asInternal().length(); }
		value_t len_sq() const noexcept { return asInternal().len_sq(); }
		bool isNaN() const noexcept { return asInternal().isNaN(); }
		bool isOutstanding() const noexcept { return asInternal().isOutstanding(); }

		auto saturation(const value_t& vMin, const value_t& vMax) const noexcept { return asInternal().saturation(vMin, vMax); }
		auto l_intp(const wrap_t& w, const value_t& r) const noexcept { return asInternal().l_intp(w, r); }
		auto absolute() const noexcept { return asInternal().absolute(); }
		value_t getMinValue() const noexcept { return asInternal().getMinValue(); }
		value_t getMaxValue() const noexcept { return asInternal().getMaxValue(); }
		void linearNormalize() noexcept { *this = linearNormalization(); }
		spec_t linearNormalization() const noexcept { return asInternal().linearNormalization(); }
		bool isZero(const value_t& th) const noexcept { return asInternal().isZero(th); }

		// ------------------- サイズ変換 -------------------
		// 余った要素はゼロで埋める
		template <int N2>
		decltype(auto) convert() const noexcept {
			return _convert<N2>(std::integral_constant<bool, (N2>size)>{});
		}
		// 大きいサイズへの変換
		template <int N2>
		auto _convert(std::true_type) const noexcept;
		// 小さいサイズへの変換
		template <int N2>
		decltype(auto) _convert(std::false_type) const noexcept;

		// 余った要素はPosの場所を任意の値、他はゼロで埋める
		template <int N2, int Pos>
		decltype(auto) convertI(const value_t& vi) const noexcept {
			return _convertI<N2,Pos>(vi, std::integral_constant<bool, (N2>size)>{});
		}
		template <int N2, int Pos>
		auto _convertI(const value_t& vi, std::true_type) const noexcept;
		template <int N2, int Pos>
		decltype(auto) _convertI(const value_t& vi, std::false_type) const noexcept;

		// -------- Luaへのエクスポート用 --------
		spec_t luaAddV(const spec_t& v) const noexcept {
			return *this + v;
		}
		spec_t luaSubV(const spec_t& v) const noexcept {
			return *this - v;
		}
		spec_t luaMulF(const float s) const noexcept {
			return *this * s;
		}
		spec_t luaMulV(const spec_t& v) const noexcept {
			return *this * v;
		}
		spec_t luaMulM(const MatT_spec<spec_t,size,size>& m) const noexcept {
			return *this * m;
		}
		spec_t luaDivF(const float s) const noexcept {
			return *this / s;
		}
		spec_t luaDivV(const spec_t& v) const noexcept {
			return *this / v;
		}
		spec_t luaInvert() const noexcept {
			return -(*this);
		}
		bool luaEqual(const spec_t& v) const noexcept {
			return *this == v;
		}
		std::string luaToString() const {
			return lubee::ToString(*this);
		}
	};
	template <class W, class D, class S>
	const int VecT<W,D,S>::size;
	template <class W, class D, class S>
	const bool VecT<W,D,S>::align;

	template <class W, class D, int N>
	struct VecT_spec : VecT<W,D, VecT_spec<W,D,N>> {
		using base_t = VecT<W,D, VecT_spec<W,D,N>>;
		using base_t::base_t;
		VecT_spec() = default;
		constexpr VecT_spec(const base_t& b) noexcept: base_t(b) {}
	};
}
#include "include/vec_d2.hpp"
#include "include/vec_d3.hpp"
#include "include/vec_d4.hpp"
namespace frea {
	template <class W, class D, class S>
	template <int N2>
	auto VecT<W,D,S>::_convert(std::true_type) const noexcept {
		return type_cn<N2>(asInternal().template convert<N2>());
	}
	template <class W, class D, class S>
	template <int N2>
	decltype(auto) VecT<W,D,S>::_convert(std::false_type) const noexcept {
		// そのままポインタ読み替え
		return reinterpret_cast<const type_cn<N2>&>(*this);
	}

	template <class W, class D, class S>
	template <int N2, int Pos>
	auto VecT<W,D,S>::_convertI(const value_t& vi, std::true_type) const noexcept {
		return type_cn<N2>(asInternal().template convertI<N2,Pos>(vi));
	}
	template <class W, class D, class S>
	template <int N2, int Pos>
	decltype(auto) VecT<W,D,S>::_convertI(const value_t&, std::false_type) const noexcept {
		// そのままポインタ読み替え
		return reinterpret_cast<const type_cn<N2>&>(*this);
	}

	template <class T, int N, bool A>
	using info_detect_t = decltype(info_detect<N,A>(std::declval<T>()));
	template <class T, int N, bool A>
	using Vec_t = info_detect_t<T,N,A>;

	template <class W, ENABLE_IF((is_wrap<W>{} || IsTuple_t<W>{}))>
	inline std::ostream& operator << (std::ostream& os, const W& w) {
		const Vec_t<typename W::value_t, W::size, true> v(w);
		return os << v;
	}
	template <class V, ENABLE_IF((is_vector<V>{}))>
	inline std::ostream& operator << (std::ostream& os, const V& v) {
		os << "Vec" << V::size << ": ";
		return v.print(os);
	}
}
namespace std {
	template <class T, int N, bool A>
	struct hash<frea::Data<T,N,A>> {
		std::size_t operator()(const frea::Data<T,N,A>& d) const noexcept {
			std::size_t ret = 0;
			for(int i=0 ; i<N ; i++)
				ret ^= *reinterpret_cast<const uint32_t*>(d.m+i);
			return ret;
		}
	};
	template <class W, class D, int N>
	struct hash<frea::VecT_spec<W,D,N>> {
		using vec_t = frea::VecT_spec<W,D,N>;
		std::size_t operator()(const vec_t& d) const noexcept {
			using base_t = typename vec_t::base_t::base_t;
			return hash<base_t>()(static_cast<const base_t&>(d));
		}
	};
}
namespace frea {
	#if SSE >= 2
		#define DEF_RV(n) \
			using RVec##n = RVec_t<float, n>;
		DEF_RV(2)
		DEF_RV(3)
		DEF_RV(4)
		#undef DEF_RV
	#endif

	#define DEF_V(n) \
		using Vec##n = Vec_t<float, n, false>; \
		using IVec##n = Vec_t<int32_t, n, false>; \
		using DVec##n = Vec_t<double, n, false>; \
		using AVec##n = Vec_t<float, n, true>;
	DEF_V(2)
	DEF_V(3)
	DEF_V(4)
	#undef DEF_V
}
