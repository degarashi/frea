#pragma once
#include "../range.hpp"
#include "vector.hpp"

namespace frea {
	namespace random {
		// --------------- matrix ---------------
		template <class V>
		constexpr V DefaultMatValue = 1e4;
		// ランダム値のデフォルト範囲
		template <class V>
		constexpr Range<V> DefaultMatRange{-DefaultMatValue<V>, DefaultMatValue<V>};
		//! 要素の値がランダムな行列
		template <class M, class RDF>
		auto GenMat(const RDF& rdf) {
			M m;
			for(int i=0 ; i<M::dim_m ; i++) {
				m.m[i] = GenVec<typename M::vec_t>(rdf);
			}
			return m;
		}
		template <class M, class RDF, class V_t=typename M::value_t>
		M GenMatRange(const RDF& rdf, const Range<V_t>& r=DefaultMatRange<V_t>) {
			return GenMat<M>([&, r](){ return rdf(r); });
		}
	}
}
