#pragma once
#include "../vector.hpp"
#include <cereal/cereal.hpp>

namespace frea {
	template <class Ar, class T, int N, bool A>
	void serialize(Ar& ar, Data<T,N,A>& d) {
		using D = Data<T,N,A>;
		std::size_t sz = D::size;
		ar(cereal::make_size_tag(sz));
		for(std::size_t i=0 ; i<sz ; i++)
			ar(d.m[i]);
	}
}
