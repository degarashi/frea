#pragma once
#include "../angle.hpp"

namespace frea {
	namespace random {
		template <class Ang, class RD>
		Ang GenAngle(RD&& rd) {
			return Ang(rd(Ang::OneRotationRange));
		}
		template <class Ang, class RD>
		Ang GenHalfAngle(RD&& rd) {
			return Ang(rd(Ang::HalfRotationRange));
		}
	}
}
