#pragma once
#include "angle.hpp"
#include "vector.hpp"
#include "../lubee/error.hpp"
#include "../quaternion.hpp"

namespace frea {
	namespace random {
		// --------------- quaternion ---------------
		//! ランダムなクォータニオン
		template <class Q, class RD>
		auto GenQuat(RD&& rd) {
			return Q::Rotation(GenVecUnit<typename Q::vec_t>(rd), GenHalfAngle<Radian<typename Q::value_t>>(rd));
		}
		template <class EQ, class RD>
		auto GenExpQuat(RD&& rd) {
			return EQ(GenQuat<QuatT<typename EQ::value_t, EQ::align>>(rd));
		}
	}
}
