#pragma once
#include <utility>
#include <tuple>
#include "meta/enable_if.hpp"
#include "meta/constant_t.hpp"

namespace frea {
	namespace seq {
		//! 2つ以上のシーケンスを結合
		/*!
			\tparam Ts	std::index_sequence群
		*/
		template <class... Ts>
		struct Concat;
		template <std::size_t... N>
		struct Concat<std::index_sequence<N...>> {
			using type = std::index_sequence<N...>;
		};
		template <std::size_t... N0, std::size_t... N1, class... Ts>
		struct Concat<std::index_sequence<N0...>, std::index_sequence<N1...>, Ts...> {
			using type = typename Concat<std::index_sequence<N0..., N1...>, Ts...>::type;
		};
		template <class... Ts>
		using Concat_t = typename Concat<Ts...>::type;

		//! 2つのstd::tupleを結合
		template <class T0, class T1>
		using TupleCat_t = decltype(std::tuple_cat(std::declval<T0>(), std::declval<T1>()));

		template <int N, int Rem>
		struct Repeat {
			using type = Concat_t<std::index_sequence<N>, typename Repeat<N,Rem-1>::type>;
		};
		template <int N>
		struct Repeat<N,0> {
			using type = std::index_sequence<>;
		};
		//! 値NをRem回繰り返したシーケンス
		template<int N, int Rem>
		using Repeat_t = typename Repeat<N,Rem>::type;

		//! シーケンスを反転
		template <class Seq>
		struct Reverse;
		template <>
		struct Reverse<std::index_sequence<>> {
			using type = std::index_sequence<>;
		};
		template <std::size_t I0, std::size_t... Idx>
		struct Reverse<std::index_sequence<I0, Idx...>> {
			using type = Concat_t<typename Reverse<std::index_sequence<Idx...>>::type, std::index_sequence<I0>>;
		};
		template <class Seq>
		using Reverse_t = typename Reverse<Seq>::type;

		template <std::size_t... N>
		auto Duplicate(SZConst<0>) -> std::index_sequence<>;
		template <std::size_t... N, int ND>
		auto Duplicate(SZConst<ND>);
		template <std::size_t... N, int ND>
		auto Duplicate(SZConst<ND>) ->
			Concat_t<
				Concat_t<
					std::index_sequence<N...>,
					std::index_sequence<N...>
				>,
				decltype(Duplicate<N...>(SZConst<ND-1>()))
			>;

		template <class Seq, std::size_t N>
		struct StripFrontIndex;
		template <std::size_t I0, std::size_t... Idx, std::size_t N>
		struct StripFrontIndex<std::index_sequence<I0, Idx...>, N> {
			using type = typename StripFrontIndex<std::index_sequence<Idx...>, N-1>::type;
		};
		template <std::size_t I0, std::size_t... Idx>
		struct StripFrontIndex<std::index_sequence<I0, Idx...>, 0> {
			using type = std::index_sequence<I0, Idx...>;
		};
		template <>
		struct StripFrontIndex<std::index_sequence<>, 0> {
			using type = std::index_sequence<>;
		};
		template <class Seq, std::size_t N>
		using StripFrontIndex_t = typename StripFrontIndex<Seq, N>::type;

		template <std::size_t I0, std::size_t... Idx>
		auto _GetFrontIndex(std::index_sequence<I0, Idx...>) -> SZConst<I0>;
		template <class Seq>
		constexpr auto GetFrontIndex = decltype(_GetFrontIndex(std::declval<Seq>())){};

		template <std::size_t... Idx>
		auto _IndexSize(std::index_sequence<Idx...>) -> SZConst<sizeof...(Idx)>;
		template <class Seq>
		constexpr auto IndexSize = decltype(_IndexSize(std::declval<Seq>())){};

		// std::tupleの先頭から指定個数分の要素を省く
		template <class Tup, std::size_t N>
		struct StripFront;
		template <class T0, class... Ts, std::size_t N>
		struct StripFront<std::tuple<T0,Ts...>, N> {
			using type = typename StripFront<std::tuple<Ts...>, N-1>::type;
		};
		template <class T0, class... Ts>
		struct StripFront<std::tuple<T0,Ts...>, 0> {
			using type = std::tuple<T0,Ts...>;
		};
		template <>
		struct StripFront<std::tuple<>, 0> {
			using type = std::tuple<>;
		};
		template <class Tup, std::size_t N>
		using StripFront_t = typename StripFront<Tup, N>::type;

		// std::tupleから指定位置の要素のみを省く
		template <std::size_t Cur, std::size_t At, class Prev, class Tail>
		struct StripAt {
			using type = typename StripAt<Cur+1, At,
							TupleCat_t<Prev, std::tuple<std::tuple_element_t<0,Tail>>>,
							StripFront_t<Tail,1>
						>::type;
		};
		template <std::size_t At, class Prev, class Tail>
		struct StripAt<At, At, Prev, Tail> {
			using type = TupleCat_t<Prev, StripFront_t<Tail, 1>>;
		};
		template <class Tup, std::size_t N>
		using StripAt_t = typename StripAt<0, N, std::tuple<>, Tup>::type;

		using BoolSeq_t = std::tuple<BConst<false>, BConst<true>>;

		template <
			template <class...> class V,
			std::size_t Rem,
			class TypeTup,
			class SeqTup
		>
		struct ExpandTypes;
		template <
			template <class...> class V,
			std::size_t Rem,
			class... Types,
			class... STypes
		>
		struct ExpandTypes<V, Rem, std::tuple<Types...>, std::tuple<STypes...>> {
			using single_type = V<std::tuple_element_t<GetFrontIndex<STypes>, Types>...>;
			using type =
				TupleCat_t<
					std::tuple<single_type>,
					typename ExpandTypes<V, Rem-1,
						std::tuple<
							Types...
						>,
						std::tuple<
							StripFrontIndex_t<STypes, 1>...
						>
					>::type
				>;
		};
		template <
			template <class...> class V,
			class... Types,
			class... STypes
		>
		struct ExpandTypes<V, 0, std::tuple<Types...>, std::tuple<STypes...>> {
			using type = std::tuple<>;
		};
		template <template <class...> class V, class TypeTup, class SeqTup>
		using ExpandTypes_t =
				typename ExpandTypes<
					V,
					IndexSize<std::tuple_element_t<0, SeqTup>>,
					TypeTup,
					SeqTup
				>::type;

		template <int From, int To>
		struct Range;
		template <int From, int To>
		struct Range {
			using type = TupleCat_t<std::tuple<IConst<From>>, typename Range<From+1,To>::type>;
		};
		template <int I>
		struct Range<I,I> {
			using type = std::tuple<>;
		};
		template <int From, int To>
		using Range_t = typename Range<From,To>::type;

		template <std::size_t From, std::size_t To>
		struct RangeIndex {
			using type = Concat_t<std::index_sequence<From>, typename RangeIndex<From+1,To>::type>;
		};
		template <std::size_t N>
		struct RangeIndex<N,N> {
			using type = std::index_sequence<>;
		};
		template <std::size_t From, std::size_t To>
		using RangeIndex_t = typename RangeIndex<From, To>::type;

		//! 特定のシーケンスを指定回ずつ、繰り返す (Tup=[0,1,2], N=3) -> [0,0,0, 1,1,1, 2,2,2]
		/*!
			\tparam Tup		対象のシーケンス (std::index_sequence)
			\tparam N		重複させる回数
		*/
		template <class Tup, std::size_t N>
		struct DuplIndex;
		template <std::size_t... Idx, std::size_t N>
		struct DuplIndex<std::index_sequence<Idx...>, N> {
			using type = Concat_t<Repeat_t<Idx, N>...>;
		};
		template <class Tup, std::size_t N>
		using DuplIndex_t = typename DuplIndex<Tup, N>::type;

		//! 特定のindex_sequenceを指定回数繰り返す
		/*!
			\tparam Tup		繰り返すシーケンス一回分(std::tuple)
			\tparam Rem		残りのリピート数
		*/
		template <class Tup, std::size_t Rem>
		struct RepeatIndex {
			using type = Concat_t<Tup, typename RepeatIndex<Tup, Rem-1>::type>;
		};
		template <class Tup>
		struct RepeatIndex<Tup, 0> {
			using type = std::index_sequence<>;
		};
		template <class Tup, std::size_t N>
		using RepeatIndex_t = typename RepeatIndex<Tup, N>::type;

		/*!
			\tparam Sum		組み合わせ総数
			\tparam Ts		組み合わせるtuple群
		*/
		template <std::size_t Sum, class... Ts>
		struct _Combi;
		// 末端処理 (T0はtuple群の末尾)
		template <std::size_t Sum, class T0>
		struct _Combi<Sum, T0> {
			constexpr static std::size_t t0_size = std::tuple_size<T0>{},					//!< この層における組み合わせ数
											size = t0_size,									//!< これを含めた下層の組み合わせ総数
											upper_size = (t0_size==0) ? 1 : Sum/t0_size;	//!< これより上層の組み合わせ総数
			// [0...t0_size-1]のシーケンスを上位層の組み合わせ回数分、繰り返す
			using index = std::tuple<RepeatIndex_t<std::make_index_sequence<t0_size>, upper_size>>;
		};
		template <std::size_t Sum, class T0, class T1, class... Ts>
		struct _Combi<Sum, T0,T1,Ts...> {
			using lower = _Combi<Sum, T1,Ts...>;
			constexpr static std::size_t lower_size = lower::size;							//!< これを含めない下層の組み合わせ総数
			using lower_index = typename lower::index;

			constexpr static std::size_t t0_size = std::tuple_size<T0>{},					//!< この層における組み合わせ数
										size = t0_size * lower_size,						//!< これを含めた下層の組み合わせ総数
										upper_size = (t0_size==0) ? 1 : Sum/size;			//!< これより上層の組み合わせ総数
			using index = TupleCat_t<
							std::tuple<
								RepeatIndex_t<
									DuplIndex_t<
										std::make_index_sequence<t0_size>,
										lower_size
									>,
									upper_size
								>
							>,
							lower_index
						>;
		};

		//! 引数の数値を全て掛けた値
		template <std::size_t... N>
		struct Mul : SZConst<1> {};
		template <std::size_t N0, std::size_t... N>
		struct Mul<N0,N...> : SZConst<N0 * Mul<N...>{}> {};

		//! 組み合わせ配列の算出
		/*!
			\tparam Tup Tupleの二次元配列 (std::tuple<std::tuple<...>, std::tuple<...>, ...>)
		*/
		template <class T>
		struct Combi;
		// 内部クラス_Combiへ、考えられる全ての組み合わせの総数と配列リストを渡す
		template <class... Ts>
		struct Combi<std::tuple<Ts...>>:
			_Combi<Mul<std::tuple_size<Ts>{}...>{}, Ts...> {};
		template <class Tup>
		using Combi_t = typename Combi<Tup>::index;

		template <
			template <class...> class V,
			class Types
		>
		using ExpandTypes_t2 = ExpandTypes_t<V, Types, Combi_t<Types>>;
	}
}
