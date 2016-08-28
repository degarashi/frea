#pragma once
#include "vector.hpp"

namespace frea {
	namespace random {
		template <class P, class RD>
		P GenPlane(RD&& rd) {
			// 中身は単位ベクトル + 原点との距離
			const auto nml = GenVecUnit<typename P::vec_t>(rd);
			const auto dist = rd();
			return P(nml, dist);
		}
	}
}
