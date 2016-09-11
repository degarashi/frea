#pragma once
#include "angle.hpp"
#include "vector.hpp"

namespace frea {
	// ------------- Angle系関数 -------------
	//! dirAを基準に反時計回りに増加する値を返す
	/*!
		\param[in] dir 値を算出したい単位ベクトル
		\param[in] dirA 基準の単位ベクトル
		\return 角度に応じた0〜4の値(一様ではない)
	*/
	template <class V, ENABLE_IF((V::size==2))>
	auto AngleValueNL(const V& dir, const V& dirA) {
		const auto d0 = dir.dot(dirA);
		if(dirA.cw(dir) <= -1e-6)
			return d0+1 + 2;
		return 2 - (d0+1);
	}
	//! 上方向を基準としたdirの角度を返す(半時計周り)
	/*!
		\return 0〜2*PI未満の角度値
	*/
	template <class V, ENABLE_IF((V::size==2))>
	auto AngleValue(const V& dir) {
		using ret_t = Radian<typename V::value_t>;
		const auto ac0 = std::acos(std::max<double>(-1.0, std::min<double>(1.0, dir.y)));
		if(dir.x >= 1e-6)
			return ret_t(ret_t::OneRotationAng - ac0);
		return ret_t(ac0);
	}
	//! ang0,ang1における-oneloop以上oneloop未満の差分角度値を計算
	template <class A>
	auto AngleLerpValueDiff(const A& ang0, const A& ang1, const A& oneloop) {
		auto diff = ang1 - ang0;
		const auto LHalf = oneloop/2;
		if(diff >= LHalf)
			diff = -oneloop + diff;
		else if(diff < -LHalf)
			diff = oneloop + diff;
		return diff;
	}
	template <class A>
	A AngleDiff(const A& ang0, const A& ang1) {
		return A(AngleLerpValueDiff(ang0.get(), ang1.get(), A::OneRotationAng));
	}

	template <class A, class Proc>
	auto AngleLerpValue(const A& ang0, const A& ang1, const Proc& proc, const A& oneloop) {
		const auto diff = AngleLerpValueDiff(ang0, ang1, oneloop);
		return ang0 + proc(diff);
	}
	//! 上方向を0度とし、反時計回り方向に指定角度回転したベクトルを返す
	inline Vec_t<float,2,false> VectorFromAngle(const RadD& ang) {
		const auto angv = ang.get();
		return {-std::sin(angv), std::cos(angv)};
	}
	//! 2つの角度値をang0 -> ang1の線形補間
	template <class T>
	T AngleLerp(const T& ang0, const T& ang1, const typename T::value_t& r) {
		const auto fn = [r](auto&& v){ return v*r; };
		return T(AngleLerpValue(ang0.get(), ang1.get(), fn, T::OneRotationAng));
	}
	//! ang0からang1へ向けてmaxDiff以下の分だけ近づける
	template <class T>
	T AngleMove(const T& ang0, const T& ang1, const T& maxDiff) {
		const auto fn = [mdiff = maxDiff.get()](auto&& v) {
			if(std::abs(v) > mdiff)
				return (v > 0) ? mdiff : -mdiff;
			return decltype(mdiff)(v);
		};
		return T(AngleLerpValue(ang0.get(), ang1.get(), fn, T::OneRotationAng));
	}
}
