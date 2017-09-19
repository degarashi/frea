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
	constexpr float PI = Pi<float>,
					SIN0 = 0,
					SIN30 = 1.f/2,
					SIN45 = 1.f/1.41421356f,
					SIN60 = 1.7320508f/2,
					SIN90 = 1.f,
					COS0 =  SIN90,
					COS30 = SIN60,
					COS45 = SIN45,
					COS60 = SIN30,
					COS90 = SIN0;
}
