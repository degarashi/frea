#pragma once
#include "../lubee/meta/enable_if.hpp"
#include "../angle.hpp"

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
		//! ランダムなベクトル (但し全ての成分の絶対値がそれぞれ基準範囲内)
		template <class V, class RD>
		auto GenVec(RD&& rd, const lubee::Range<typename V::value_t>& rTh) {
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
	}
}
