#pragma once
#include "../affine_parts.hpp"
#include "vector.hpp"

namespace frea {
	template <class Ar, class T>
	void serialize(Ar& ar, AffineParts<T>& af) {
		ar(
			cereal::make_nvp("offset", af.offset),
			cereal::make_nvp("scale", af.scale),
			cereal::make_nvp("rotation", af.rotation)
		);
	}
}
