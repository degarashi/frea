#pragma once
#include "../range.hpp"

namespace frea {
	namespace random {
		template <class T, class RDF>
		Range<T> GenRange(const RDF& rdf) {
			T rmin = rdf(),
			 rmax = rdf();
			if(rmin > rmax)
				std::swap(rmin, rmax);
			return {rmin, rmax};
		}
	}
}
