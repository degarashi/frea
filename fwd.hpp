#pragma once

namespace frea {
	// (affine_parts.hpp)
	template <class T>
	struct AffineParts;

	// (angle.hpp)
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

	// (vector.hpp)
	template <class T, int N, bool A>
	struct Data;
	template <class W, class D, int N>
	struct VecT_spec;

	// (matrix.hpp)
	template <class V, int M>
	class DataM;
	template <class V, int M, int N>
	struct MatT_spec;

	// (quaternion.hpp)
	template <class T, bool A>
	struct QuatT;
	using Quat = QuatT<float, false>;
	using AQuat = QuatT<float, true>;
	using DQuat = QuatT<double, false>;
	using ADQuat = QuatT<double, true>;

	// (expquat.hpp)
	template <class T, bool A>
	struct ExpQuatT;
	using ExpQuat = ExpQuatT<float, false>;
	using AExpQuat = ExpQuatT<float, true>;
	using DExpQuat = ExpQuatT<double, false>;
	using ADExpQuat = ExpQuatT<double, true>;

	// (plane.hpp)
	template <class T, bool A>
	struct PlaneT;
	using Plane = PlaneT<float, false>;
	using APlane = PlaneT<float, true>;
	using DPlane = PlaneT<double, false>;
	using ADPlane = PlaneT<double, true>;
}
