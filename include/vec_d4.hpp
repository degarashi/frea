//! 4次元ベクトル特有の関数など定義
//! (vector.hppからインクルード)
#pragma once

namespace frea {
	template <class W, class D>
	struct VecT_spec<W, D, 4> : VecT<W,D, VecT_spec<W,D,4>> {
		using base_t = VecT<W,D, VecT_spec<W,D,4>>;
		using base_t::base_t;
		using wrap_t = typename base_t::wrap_t;

		using this_t = VecT_spec;
		typename base_t::template type_cn<3>
			asVec3Coord() const { return this->asInternal().asVec3Coord(); }
	};
}
