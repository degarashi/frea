#pragma once
#include "../yawpitchdist.hpp"
#include <cereal/access.hpp>

namespace frea {
	template <class Ar, class T>
	void serialize(Ar& ar, YawPitchDist<T>& y) {
		ar(
			cereal::make_nvp("yaw", y.yaw),
			cereal::make_nvp("pitch", y.pitch),
			cereal::make_nvp("distance", y.distance)
		);
	}
}
