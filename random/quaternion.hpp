#pragma once
#include "../error.hpp"
#include "../angle.hpp"
#include "../quaternion.hpp"
#include "vector.hpp"

namespace frea {
	namespace random {
		// --------------- quaternion ---------------
		//! ランダムなクォータニオン
		template <class Q, class RD>
		auto GenQuat(RD&& rd) {
			return Q::Rotation(GenDir<typename Q::vec_t>(rd), Radian<typename Q::value_t>(rd));
		}
		template <class EQ, class RD>
		auto GenExpQuat(RD&& rd) {
			return EQ(GenQuat<QuatT<typename EQ::value_t, EQ::align>>(rd));
		}
	}
}
