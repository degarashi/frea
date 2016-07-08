#pragma once
#include <utility>

namespace frea {
	template <int N>
	using IConst = std::integral_constant<int, N>;
	template <bool A>
	using BConst = std::integral_constant<bool, A>;
	template <std::size_t N>
	using SZConst = std::integral_constant<std::size_t, N>;
	template <class T, T val0, T val1>
	using integral_pair = std::pair<std::integral_constant<T, val0>,
									std::integral_constant<T, val1>>;
	template <int N0, int N1>
	using IPair = integral_pair<int, N0,N1>;
}
