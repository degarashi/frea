#pragma once
#include "detect_type.hpp"
#include "lubee/meta/enable_if.hpp"

namespace frea {
	template <class OItr, class Itr0, class Itr1, class T>
	void Lerp(
		OItr dst,
		Itr0 itr0,
		const Itr0 itr0E,
		Itr1 itr1,
		const T& t
	) {
		while(itr0 != itr0E) {
			*dst = *itr0 + (*itr1 - *itr0) * t;
			++dst; ++itr0; ++itr1;
		}
	}
	template <class V, class T, ENABLE_IF(!is_quaternion<V>{})>
	auto Lerp(const V& v0, const V& v1, const T& t) noexcept {
		return v0 + (v1-v0)*t;
	}
	template <class V, class T, ENABLE_IF(is_quaternion<V>{})>
	auto Lerp(const V& v0, const V& v1, const T& t) noexcept {
		return v0.slerp(v1, t);
	}
}
