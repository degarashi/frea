#pragma once
#include <type_traits>
#include <limits>
#include <cmath>

namespace frea {
	//! コンパイル時数値計算 & 比較
	template <int M, int N>
	struct Arithmetic {
		enum {
			add = M+N,
			sub = M-N,
			less = (M>N) ? N : M,
			great = (M>N) ? M : N,
			less_eq = (M<=N) ? 1 : 0,
			great_eq = (M>=N) ? 1 : 0,
			lesser = (M<N) ? 1 : 0,
			greater = (M>N) ? 1 : 0,
			equal = M==N ? 1 : 0
		};
	};
	//! bool -> std::true_type or std::false_type
	template <int V>
	using TFCond = typename std::conditional<V!=0, std::true_type, std::false_type>::type;
	//! 2つの定数の演算結果をstd::true_type か std::false_typeで返す
	template <int M, int N>
	struct ArithmeticT {
		using t_and = TFCond<(M&N)>;
		using t_or = TFCond<(M|N)>;
		using t_xor = TFCond<(M^N)>;
		using t_nand = TFCond<((M&N) ^ 0x01)>;
		using less = TFCond<(M<N)>;
		using great = TFCond<(M>N)>;
		using equal = TFCond<M==N>;
		using not_equal = TFCond<M!=N>;
		using less_eq = TFCond<(M<=N)>;
		using great_eq = TFCond<(M>=N)>;
	};
	//! 2つのintegral_constant<bool>を論理演算
	template <class T0, class T1>
	struct Logistic {
		constexpr static int I0 = std::is_same<T0, std::true_type>::value,
							I1 = std::is_same<T1, std::true_type>::value;
		using t_and = TFCond<I0&I1>;
		using t_or = TFCond<I0|I1>;
		using t_xor = TFCond<I0^I1>;
		using t_nand = TFCond<(I0&I1) ^ 0x01>;
	};
	//! 最大値を取得
	template <int... N>
	struct IMax {
		constexpr static int result = 0;
	};
	template <int N0, int... N>
	struct IMax<N0, N...> {
		constexpr static int result = N0;
	};
	template <int N0, int N1, int... N>
	struct IMax<N0, N1, N...> {
		constexpr static int result = Arithmetic<N0, IMax<N1, N...>::result>::great;
	};
}
