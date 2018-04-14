#pragma once
#include "detect_type.hpp"
#include "lubee/meta/enable_if.hpp"

namespace frea {
	template <class V, class T, ENABLE_IF(!is_quaternion<V>{})>
	auto Lerp(const V& v0, const V& v1, const T& t) noexcept {
		return v0 + (v1-v0)*t;
	}
	template <class V, class T, ENABLE_IF(is_quaternion<V>{})>
	auto Lerp(const V& v0, const V& v1, const T& t) noexcept {
		return v0.slerp(v1, t);
	}
}
