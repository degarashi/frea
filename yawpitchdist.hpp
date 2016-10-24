#pragma once
#include "quaternion.hpp"

namespace frea {
	template <class T>
	struct YawPitchDist : lubee::op::Ne<YawPitchDist<T>> {
		using Rad = Radian<T>;
		using Vec3 = Vec_t<T,3,true>;
		using Quat = QuatT<T,true>;
		Rad		yaw, pitch;
		T		distance;

		template <class Ar, class T2>
		friend void serialize(Ar&, YawPitchDist<T2>&);

		bool operator == (const YawPitchDist& ypd) const noexcept {
			return yaw==ypd.yaw &&
					pitch==ypd.pitch &&
					distance==ypd.distance;
		}
		//! 方向ベクトルをYaw,Pitch,Distanceに分解
		template <class V, ENABLE_IF((is_vector<V>{} || is_wrap<V>{}))>
		static auto FromPos(const V& pos) noexcept {
			YawPitchDist<T> ypd;
			// Distance
			Vec3 v(pos);
			ypd.distance = v.length();
			v /= ypd.distance;

			// Yaw
			Vec3 xzvec(v.x, 0, v.z);
			if(xzvec.len_sq() < lubee::ThresholdF<T>(0.4))
				ypd.yaw.set(0);
			else {
				xzvec.normalize();
				T ac = std::acos(lubee::Saturate<T>(xzvec.z, 1));
				if(xzvec.x < 0)
					ac = 2*Pi<T> - ac;
				ypd.yaw.set(ac);
			}

			// Pitch
			constexpr T h = AngleInfo<Radian_t>::one_rotation<T> / 2;
			ypd.pitch.set(lubee::Saturate(std::asin(lubee::Saturate<T>(v.y, 1)), -h, h));
			return ypd;
		}
		//! YawPitchDistの位置から座標原点を見る姿勢
		auto toOffsetRot() const {
			struct {
				Vec3	pos;
				Quat	rot;
			} ret;
			// Z軸をYaw/Pitch/Roll傾けた方向に対してDist距離進んだ場所がカメラの位置
			// カメラの方向は変換済みZ軸と逆
			const Quat q = Quat::RotationYPR(yaw, pitch, RadF(0));
			const Vec3 z = q.getDir();
			ret.pos = z*distance;
			const Vec3 vd = -(z*distance).normalization();
			ret.rot = Quat::LookAt(vd, Vec3(0,1,0));
			return ret;
		}
	};
	template <class T>
	inline std::ostream& operator << (std::ostream& os, const YawPitchDist<T>& ypd) {
		return os << "YPD: yaw=" << ypd.yaw << " pitch=" << ypd.pitch << " dist=" << ypd.distance;
	}
}
