#pragma once
#include "enable_if.hpp"
#include <type_traits>
#include <utility>

namespace frea {
	//! コンパイル時定数で数値のN*10乗を計算
	template <class T, int N, ENABLE_IF((N==0))>
	constexpr T ConstantPow10(T value=1, IConst<N>* =nullptr) {
		return value;
	}
	template <class T, int N, ENABLE_IF((N<0))>
	constexpr T ConstantPow10(T value=1, IConst<N>* =nullptr) {
		return ConstantPow10(value/10, (IConst<N+1>*)nullptr);
	}
	template <class T, int N, ENABLE_IF((N>0))>
	constexpr T ConstantPow10(T value=1, IConst<N>* =nullptr) {
		return ConstantPow10(value*10, (IConst<N-1>*)nullptr);
	}
}
