#pragma once
#include "quaternion.hpp"

namespace frea {
	template <class T>
	struct AffineParts : lubee::op::Ne<AffineParts<T>> {
		using Vec = Vec_t<T,3,true>;
		using Quat = QuatT<T,true>;
		Vec		offset,
				scale;
		Quat	rotation;

		template <class Ar, class T2>
		friend void serialize(Ar&, AffineParts<T2>&);

		bool operator == (const AffineParts& ap) const noexcept {
			return offset == ap.offset &&
				scale == ap.scale &&
				rotation == ap.rotation;
		}
		//! アフィン成分分解
		template <class M, ENABLE_IF((is_matrix<M>{} || is_wrapM<M>{}) && (M::dim_m==4 && M::dim_n>=3))>
		static AffineParts<T> Decomp(const M& m) noexcept {
			AffineParts<T> ap;
			Mat_t<T, 3, 3, true> tm;
			// オフセットは4行目をそのまま抽出
			ap.offset = m.template getRow<3>().template convert<3>();
			for(int i=0 ; i<3 ; i++) {
				const auto r = m[i].template convert<3>();
				ap.scale[i] = r.length();
				tm[i] = r / ap.scale[i];
			}
			// 掌性チェックして、左手系でなければX軸を反転
			const auto cv = tm.template getRow<1>().cross(tm.template getRow<2>());
			if(tm.template getRow<0>().dot(cv) <= 0) {
				tm.template getRow<0>() *= -1;
				ap.scale.x *= -1;
			}
			ap.rotation = Quat::FromAxis(tm.template getRow<0>(), tm.template getRow<1>(), tm.template getRow<2>());
			return ap;
		}
	};
	template <class T>
	inline std::ostream& operator << (std::ostream& os, const AffineParts<T>& ap) {
		os << "AffineParts: offset=";
		ap.offset.print(os);
		os << ", scale=";
		ap.scale.print(os);
		os << ", rotation=";
		ap.rotation.print(os);
		return os;
	}
}
