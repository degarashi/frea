#pragma once

namespace frea {
	struct Axis {
		enum e {
			X,
			Y,
			Z,
			_Num
		};
	};
	template <class T>
	T Square(const T& t) {
		return t*t;
	}
	template <class T>
	T Cubic(const T& t) {
		return t*t*t;
	}
	template <class T>
	constexpr T Pi = T(3.1415926535);	// std::atan(1.0f)*4
}
