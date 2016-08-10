//! 2次元ベクトル特有の関数など定義
//! (vector.hppからインクルード)
#pragma once

namespace frea {
	template <class W, class D>
	struct VecT_spec<W, D, 2> : VecT<W,D, VecT_spec<W,D,2>> {
		using base_t = VecT<W,D, VecT_spec<W,D,2>>;
		using base_t::base_t;
		using wrap_t = typename base_t::wrap_t;
		using this_t = VecT_spec;
		using value_t = typename base_t::value_t;

		value_t ccw(const wrap_t& v) const { return this->asInternal().ccw(v); }
		value_t cw(const wrap_t& v) const { return this->asInternal().cw(v); }
		static value_t Ccw(const base_t& v0, const base_t& v1, const base_t& v2) {
			return (v0-v1).ccw(v2-v1);
		}
		static value_t Cw(const base_t& v0, const base_t& v1, const base_t& v2) {
			return (v0-v1).cw(v2-v1);
		}
	};
}
