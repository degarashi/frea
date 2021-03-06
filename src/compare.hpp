#pragma once
#include "lubee/src/meta/enable_if.hpp"
#include "lubee/src/meta/constant.hpp"
#include <limits>
#include <cmath>

namespace frea {
	// 値比較(1次元: m[N])
	template <class T, class CMP, int N>
	bool _EqFunc(const T& v0, const T& v1, CMP&& cmp, lubee::integral_pair<int,0,N>) {
		for(int i=0 ; i<N ; i++) {
			if(!cmp(v0.m[i], v1.m[i]))
				return false;
		}
		return true;
	}
	// 値比較(2次元: ma[M][N])
	template <class T, class CMP, int M, int N>
	bool _EqFunc(const T& v0, const T& v1, CMP&& cmp, lubee::integral_pair<int,M,N>) {
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
	bool _EqFunc(const T& v0, const T& v1, CMP&& cmp, lubee::integral_pair<int,0,0>) {
		return cmp(v0, v1);
	}

	//! 絶対値の誤差による等値判定
	/*! \param[in] val value to check
		\param[in] vExcept target value
		\param[in] vEps value threshold */
	template <class T, class T2>
	bool EqAbs(const T& val, const T& vExcept, T2 vEps = std::numeric_limits<T>::epsilon()) {
		auto fnCmp = [vEps](const auto& val, const auto& except){ return std::abs(except-val) <= vEps; };
		return fnCmp(val, vExcept);
	}
	template <class T, class... Ts>
	bool EqAbsT(const std::tuple<Ts...>& /*tup0*/, const std::tuple<Ts...>& /*tup1*/, const T& /*epsilon*/, lubee::IConst<-1>*) {
		return true;
	}
	//! std::tuple全部の要素に対してEqAbsを呼ぶ
	template <class T, int N, class... Ts, ENABLE_IF((N>=0))>
	bool EqAbsT(const std::tuple<Ts...>& tup0, const std::tuple<Ts...>& tup1, const T& epsilon, lubee::IConst<N>*) {
		return EqAbs(std::get<N>(tup0), std::get<N>(tup1), epsilon)
				&& EqAbsT(tup0, tup1, epsilon, (lubee::IConst<N-1>*)nullptr);
	}
	template <class... Ts, class T>
	bool EqAbs(const std::tuple<Ts...>& tup0, const std::tuple<Ts...>& tup1, const T& epsilon) {
		return EqAbsT<T>(tup0, tup1, epsilon, (lubee::IConst<sizeof...(Ts)-1>*)nullptr);
	}

	//! std::tupleの要素ごとの距離(EqAbs)比較
	template <class T, int NPow>
	struct TupleNear {
		template <class P>
		bool operator()(const P& t0, const P& t1) const {
			return EqAbs(t0, t1, lubee::ConstantPow<NPow>(T(10)));
		}
	};
}
