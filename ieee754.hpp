#pragma once
#include "meta/enable_if.hpp"
#include <limits>

namespace frea {
	//! 指数部や仮数部に使用されるビット数を定義
	template <class T>
	struct IEEE754;
	struct Half;

	template <>
	struct IEEE754<Half> {
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

	//! 浮動小数点数の精度に応じた許容誤差を算出 (第二引数は無視)
	/*!
		\param fr 仮数部を切り捨てる割合(0〜1)
	*/
	template <class T,
		 ENABLE_IF((
			std::is_same<T,Half>{} ||
			std::is_floating_point<T>{}
	))>
	constexpr T Threshold(const double fr, const T&) noexcept {
		constexpr auto FB = IEEE754<T>::FracBits;
		const auto r = int64_t(1) << static_cast<int>(FB*fr);
		return std::numeric_limits<T>::epsilon()*r;
	}
	//! 浮動小数点数の精度に応じた許容誤差を算出 (整数版 = 第二引数がそのまま返る)
	template <class T, ENABLE_IF(std::is_integral<T>{})>
	constexpr T Threshold(const double, const T& ia) noexcept { return ia; }

	template <class T,
		 ENABLE_IF((
			std::is_same<T,Half>{} ||
			std::is_floating_point<T>{}
	))>
	constexpr T ThresholdF(const double fr) noexcept {
		return Threshold<T>(fr, 0);
	}
}
