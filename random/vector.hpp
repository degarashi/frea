#pragma once
#include "../range.hpp"

namespace frea {
	namespace random {
		//! 要素の値がランダムなベクトル
		/*!
			\tparam V
			\tparam RDF
		*/
		template <class V, class RD>
		auto GenVec(RD&& rd) {
			V ret;
			for(auto& v : ret.m)
				v = rd();
			return ret;
		}
		//! 条件を満たすランダムベクトルを生成
		template <class V, class RD, class Chk>
		auto GenVecCnd(RD&& rd, const Chk& chk) {
			V v;
			do {
				v = GenVec<V>(rd);
			} while(!chk(v));
			return v;
		}
		//! ランダムなベクトル（但し長さがゼロではない）
		template <class V, class RD>
		auto GenVecLen(RD&& rd, const typename V::value_t& th) {
			return GenVecCnd<V>(rd, [th](const auto& v){ return v.length() >= th; });
		}
		//! ランダムなベクトル (但し全ての成分の絶対値がそれぞれ基準範囲内)
		template <class V, class RD>
		auto GenVecAbs( RD&& rd, const Range<typename V::value_t>& rTh) {
			return GenVecCnd<V>(rd,
					[rTh](const auto& v){
						for(const auto& vm : v.m) {
							if(!rTh.hit(vm))
								return false;
						}
						return true;
					});
		}
		//! ランダムな方向ベクトル
		template <class V, class RD>
		auto GenDir(RD&& rd) {
			using v_t = typename V::value_t;
			return
				GenVecLen<V>(
					rd,
					{-1e2, 1e2},
					{-v_t(1), v_t(1)}
				).normalization();
		}
	}
}
