#pragma once
#include "../angle.hpp"

namespace frea {
	namespace random {
		template <class Ang, class RD>
		Ang GenAngle(RD&& rd) {
			constexpr auto OR = Ang::OneRotationRange;
			return Ang(rd(OR));
		}
		template <class Ang, class RD>
		Ang GenHalfAngle(RD&& rd) {
			constexpr auto HR = Ang::HalfRotationRange;
			return Ang(rd(HR));
		}
	}
}
