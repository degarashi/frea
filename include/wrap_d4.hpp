//! 4次元ベクトル特有の関数など定義
//! (vector.hppからインクルード)
#pragma once

namespace frea {
	template <class R>
	struct wrap_spec<R,4> : wrap<R,4,wrap_spec<R,4>> {
		using base_t = wrap<R,4,wrap_spec<R,4>>;
		using base_t::base_t;

		auto asVec3Coord() const {
			auto tmp = *this;
			tmp.template makeEquality<3>();
			return this->template convert<3>() * base_t::I::Reciprocal(tmp);
		}
	};
}
