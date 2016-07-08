//! 3次元ベクトル特有の関数など定義
//! (vector.hppからインクルード)
#pragma once

namespace frea {
	template <class W, class D>
	struct VecT_spec<W, D, 3> : VecT<W,D, VecT_spec<W,D,3>> {
		using base_t = VecT<W,D, VecT_spec<W,D,3>>;
		using base_t::base_t;
		using wrap_t = typename base_t::wrap_t;

		using this_t = VecT_spec;
		this_t cross(const wrap_t& v) const { return this->::asInternal().cross(v); }
		this_t verticalVector() const { return this->asInternal().verticalVector(); }
	};
}
