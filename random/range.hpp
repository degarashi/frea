#pragma once
#include "../range.hpp"

namespace frea {
	namespace random {
		template <class T, class RD>
		Range<T> GenRange(RD&& rd) {
			auto	rmin = rd(),
					rmax = rd();
			if(rmin > rmax)
				std::swap(rmin, rmax);
			return {rmin, rmax};
		}
	}
}
