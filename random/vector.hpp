#pragma once
#include "../range.hpp"

namespace frea {
	namespace random {
		//! 要素の値がランダムなベクトル
		/*!
			\tparam V
			\tparam RDF
		*/
		template <class V, class RDF>
		auto GenVec(const RDF& rdf) {
			V ret;
			for(auto& v : ret.m)
				v = rdf();
			return ret;
		}
		//! ランダム値のデフォルト範囲値
		template <class V>
		constexpr V DefaultVecValue = 1e4;
		//! ランダム値のデフォルト範囲
		template <class V>
		constexpr Range<V> DefaultVecRange{-DefaultVecValue<V>, DefaultVecValue<V>};
		template <class V, class RDF, class V_t=typename V::value_t>
		auto GenVecRange(const RDF& rdf, const Range<V_t>& r=DefaultVecRange<V_t>) {
			return GenVec<V>([&rdf, r](){ return rdf(r); });
		}
		//! 条件を満たすランダムベクトルを生成
		template <class V, class RDF, class Chk, class V_t=typename V::value_t>
		auto GenVecCnd(const RDF& rdf, const Chk& chk, const Range<V_t>& r=DefaultVecRange<V_t>) {
			V v;
			do {
				v = GenVecRange<V>(rdf, r);
			} while(!chk(v));
			return v;
		}
		//! ランダムなベクトル（但し長さがゼロではない）
		template <class V, class RDF, class V_t=typename V::value_t>
		auto GenVecLen(const RDF& rdf,
			const V_t& th,
			const Range<V_t>& r=DefaultVecRange<V_t>)
		{
			return GenVecCnd<V>(rdf, [th](const auto& v){ return v.length() >= th; }, r);
		}
		//! ランダムなベクトル (但し全ての成分の絶対値がそれぞれ基準範囲内)
		template <class V, class RDF, class V_t=typename V::value_t>
		auto GenVecAbs(
			const RDF& rdf,
			const Range<V_t>& rTh,
			const Range<V_t>& r=DefaultVecRange<V_t>
		) {
			return GenVecCnd<V>(rdf,
					[rTh](const auto& v){
						for(const auto& vm : v.m) {
							if(!rTh.hit(vm))
								return false;
						}
						return true;
					},
					r);
		}
		//! ランダムな方向ベクトル
		template <class V, class RDF>
		auto GenDir(const RDF& rdf) {
			using v_t = typename V::value_t;
			return
				GenVecLen<V>(
					rdf,
					DefaultVecValue<v_t>,
					{-v_t(1), v_t(1)}
				).normalization();
		}
	}
}
