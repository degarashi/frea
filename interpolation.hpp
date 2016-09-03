#pragma once

namespace frea {
	template <class T>
	T Lerp(const T& t0, const T& t1, const T& t) noexcept {
		return t0 + (t1-t0)*t;
	}
}
