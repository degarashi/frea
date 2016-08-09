#pragma once
#include "../range.hpp"
#include "vector.hpp"

namespace frea {
	namespace random {
		// --------------- matrix ---------------
		//! 要素の値がランダムな行列
		template <class M, class RD>
		auto GenMat(RD&& rd) {
			M m;
			for(int i=0 ; i<M::dim_m ; i++) {
				m.m[i] = GenVec<typename M::vec_t>(rd);
			}
			return m;
		}
	}
}
