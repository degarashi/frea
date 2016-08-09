#pragma once
#include "vector.hpp"

namespace frea {
	namespace random {
		// --------------- quaternion ---------------
		//! ランダムなクォータニオン
		template <class Q, class RDF>
		auto GenQuat(const RDF& rdf) {
			return Q::Rotation(GenDir<typename Q::vec_t>(rdf), Radian<typename Q::value_t>(rdf));
		}
		template <class EQ, class RDF>
		auto GenExpQuat(const RDF& rdf) {
			return EQ(GenQuat<QuatT<typename EQ::value_t, EQ::align>>(rdf));
		}
	}
}
