#pragma once

namespace frea {
	struct Degree_t;
	struct Radian_t;

	template <class TAG, class V>
	class Angle;
	template <class T>
	using Degree = Angle<Degree_t, T>;
	using DegF = Degree<float>;
	using DegD = Degree<double>;
	template <class T>
	using Radian = Angle<Radian_t, T>;
	using RadF = Radian<float>;
	using RadD = Radian<double>;
}
