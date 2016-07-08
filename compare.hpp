#pragma once
#include "meta/enable_if.hpp"
#include "meta/constant.hpp"
#include "meta/size.hpp"
#include <limits>
#include <cmath>

namespace frea {
	// 値比較(1次元: m[N])
	template <class T, class CMP, int N>
	bool _EqFunc(const T& v0, const T& v1, CMP&& cmp, integral_pair<int,0,N>) {
		for(int i=0 ; i<N ; i++) {
			if(!cmp(v0.m[i], v1.m[i]))
				return false;
		}
		return true;
	}
	// 値比較(2次元: ma[M][N])
	template <class T, class CMP, int M, int N>
	bool _EqFunc(const T& v0, const T& v1, CMP&& cmp, integral_pair<int,M,N>) {
		for(int i=0 ; i<M ; i++) {
			for(int j=0 ; j<N ; j++) {
				if(!cmp(v0.ma[i][j], v1.ma[i][j]))
					return false;
			}
		}
		return true;
	}
	// 値比較(単一値)
	template <class T, class CMP>
	bool _EqFunc(const T& v0, const T& v1, CMP&& cmp, integral_pair<int,0,0>) {
		return cmp(v0, v1);
	}

	//! 絶対値の誤差による等値判定
	/*! \param[in] val value to check
		\param[in] vExcept target value
		\param[in] vEps value threshold */
	template <class T, class T2>
	bool EqAbs(const T& val, const T& vExcept, T2 vEps = std::numeric_limits<T>::epsilon()) {
		auto fnCmp = [vEps](const auto& val, const auto& except){ return std::fabs(except-val) <= vEps; };
		return _EqFunc(val, vExcept, fnCmp, decltype(GetWidthHeightT<T>())());
	}
	template <class T, class... Ts>
	bool EqAbsT(const std::tuple<Ts...>& /*tup0*/, const std::tuple<Ts...>& /*tup1*/, const T& /*epsilon*/, IConst<-1>*) {
		return true;
	}
	//! std::tuple全部の要素に対してEqAbsを呼ぶ
	template <class T, int N, class... Ts, ENABLE_IF((N>=0))>
	bool EqAbsT(const std::tuple<Ts...>& tup0, const std::tuple<Ts...>& tup1, const T& epsilon, IConst<N>*) {
		return EqAbs(std::get<N>(tup0), std::get<N>(tup1), epsilon)
				&& EqAbsT(tup0, tup1, epsilon, (IConst<N-1>*)nullptr);
	}
	template <class... Ts, class T>
	bool EqAbs(const std::tuple<Ts...>& tup0, const std::tuple<Ts...>& tup1, const T& epsilon) {
		return EqAbsT<T>(tup0, tup1, epsilon, (IConst<sizeof...(Ts)-1>*)nullptr);
	}

	//! 浮動少数点数の値がNaNになっているか
	template <class T, ENABLE_IF((std::is_floating_point<T>{}))>
	bool IsNaN(const T& val) {
		return !(val>=T(0)) && !(val<T(0));
	}
	//! 浮動少数点数の値がNaN又は無限大になっているか
	template <class T, ENABLE_IF((std::is_floating_point<T>{}))>
	bool IsOutstanding(const T& val) {
		auto valA = std::fabs(val);
		return valA==std::numeric_limits<T>::infinity() || IsNaN(valA);
	}
	//! 値飽和
	template <class T>
	T Saturate(const T& val, const T& minV, const T& maxV) {
		if(val > maxV)
			return maxV;
		if(val < minV)
			return minV;
		return val;
	}
	template <class T>
	T Saturate(const T& val, const T& range) {
		return Saturate(val, -range, range);
	}
	//! 値の範囲判定
	template <class T>
	bool IsInRange(const T& val, const T& vMin, const T& vMax) {
		return val>=vMin && val<=vMax;
	}
	template <class T>
	bool IsInRange(const T& val, const T& vMin, const T& vMax, const T& vEps) {
		return IsInRange(val, vMin-vEps, vMax+vEps);
	}
	//! std::tupleの要素ごとの距離(EqAbs)比較
	template <class T, int NPow>
	struct TupleNear {
		template <class P>
		bool operator()(const P& t0, const P& t1) const {
			return EqAbs(t0, t1, frea::ConstantPow10<T,NPow>());
		}
	};
}
