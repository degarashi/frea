#pragma once
#include <type_traits>
#include "constant_t.hpp"
#include "enable_if.hpp"

namespace frea {
	template <class T, int N>
	constexpr T _BitPaste(const T& val, IConst<N>, IConst<N>) {
		return val;
	}
	template <class T, int W, int N>
	constexpr T _BitPaste(const T& val, IConst<W>, IConst<N>) {
		return _BitPaste(val|(val>>W), IConst<W*2>(), IConst<N>());
	}
	//! 値が負数なら全てのビットが立った数、そうでなければゼロを返す
	template <class T, ENABLE_IF((std::is_integral<T>{}))>
	constexpr T SignMask(const T& val) {
		constexpr int S = sizeof(T)*8;
		constexpr T Mask(T(1)<<(S-1));
		return _BitPaste(val&Mask, IConst<1>(), IConst<sizeof(T)*8>());
	}
}
