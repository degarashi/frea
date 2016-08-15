#pragma once
#include <cstdint>

namespace frea {
	//! 指数部や仮数部に使用されるビット数を定義
	template <class T>
	struct IEEE754;

	// Half
	template <>
	struct IEEE754<short> {
		using Integral_t = int16_t;
		constexpr static int ExpBits = 5,
							FracBits = 10,
							ExpZero = 15;
	};
	template <>
	struct IEEE754<float> {
		using Integral_t = int32_t;
		constexpr static int ExpBits = 8,
							FracBits = 23,
							ExpZero = 127;
	};
	template <>
	struct IEEE754<double> {
		using Integral_t = int64_t;
		constexpr static int ExpBits = 11,
							FracBits = 52,
							ExpZero = 1023;
	};
}
