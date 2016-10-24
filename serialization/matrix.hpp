#pragma once
#include "../matrix.hpp"
#include "vector.hpp"
#include <cereal/cereal.hpp>

namespace frea {
	template <class Ar, class V, int M>
	void serialize(Ar& ar, DataM<V,M>& d) {
		using Dm = DataM<V,M>;
		std::size_t sz = Dm::dim_m;
		ar(cereal::make_size_tag(sz));
		for(std::size_t i=0 ; i<sz ; i++)
			ar(d.m[i]);
	}
}
