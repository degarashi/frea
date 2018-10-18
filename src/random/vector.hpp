#pragma once
#include "lubee/src/meta/enable_if.hpp"
#include "../angle.hpp"
#include <vector>

namespace frea {
	namespace random {
		//! 要素の値がランダムなベクトル
		/*!
			\tparam V
			\tparam RDF
		*/
		template <class V, class RD>
		V GenVec(RD&& rd) {
			V ret;
			for(auto& v : ret.m)
				v = rd();
			return ret;
		}
		//! ランダムなベクトル (但し全ての成分の絶対値がそれぞれ基準範囲内)
		template <class V, class RD>
		V GenVec(RD&& rd, const lubee::Range<typename V::value_t>& rTh) {
			V ret;
			for(auto& v : ret.m)
				v = rd(rTh);
			return ret;
		}
		//! ランダムな単位ベクトル(2次元)
		template <class V, class RD, ENABLE_IF(V::size==2)>
		V GenVecUnit(RD&& rd) {
			using value_t = typename V::value_t;
			static_assert(std::is_floating_point<value_t>{}, "");
			constexpr auto OR = Radian<value_t>::OneRotationRange;
			const auto ang = rd(OR);
			return {std::sin(ang), std::cos(ang)};
		}
		//! ランダムな単位ベクトル(3次元)
		template <class V, class RD, ENABLE_IF(V::size==3)>
		V GenVecUnit(RD&& rd) {
			using value_t = typename V::value_t;
			static_assert(std::is_floating_point<value_t>{}, "");
			constexpr auto OR = Radian<value_t>::OneRotationRange;
			const auto yaw = rd(OR),
						pitch = std::asin(rd({-1,1}));
			return {
				std::cos(yaw) * std::cos(pitch),
				std::sin(pitch),
				std::sin(yaw) * std::cos(pitch)
			};
		}
		namespace detail {
			template <class V, class Gen, class Chk>
			auto MakeVec(const int n, Gen&& gen, Chk&& chk) {
				D_Assert(n>0,  "n shoud be greater than 0");
				std::vector<V> ret(n);
				ret[0] = gen();
				for(int i=1 ; i<n ; i++) {
					for(;;) {
						ret[i] = gen();
						bool b = true;
						for(int j=0 ; j<i ; j++) {
							if(!chk(ret[j], ret[i])) {
								b = false;
								break;
							}
						}
						if(b)
							break;
					}
				}
				return ret;
			}
		}
		//! 方向が重ならない単位ベクトルを任意の数、生成
		template <class V, class RDF>
		std::vector<V> GenVecUnitN(RDF&& rdf, const int n, const typename V::value_t& th) {
			return detail::MakeVec<V>(
				n,
				[&rdf](){ return GenVecUnit<V>(rdf); },
				[th](const auto& v0, const auto& v1){ return v0.dot(v1) < th; }
			);
		}
		//! 位置が重ならない座標ベクトルを任意の数、生成
		template <class V, class RDF>
		std::vector<V> GenVecN(RDF&& rdf, const int n, const typename V::value_t& th) {
			return detail::MakeVec<V>(
				n,
				[&rdf](){ return GenVec<V>(rdf); },
				[th=th*th](const auto& v0, const auto& v1){ return v0.dist_sq(v1) >= th; }
			);
		}
	}
}
