//! 3次元ベクトル特有の関数など定義
//! (vector.hppからインクルード)
#pragma once
#include <tuple>

namespace frea {
	template <class R>
	struct wrap_spec<R,3> : wrap<R,3,wrap_spec<R,3>> {
		using base_t = wrap<R,3,wrap_spec<R,3>>;
		using base_t::base_t;
		using this_t = wrap_spec;
		using value_t = typename base_t::value_t;

		using base_t::operator *;
		template <class Q,
				 ENABLE_IF(is_quaternion<Q>{})>
		auto operator * (const Q& q) const {
			QuatT<typename Q::value_t, true> q0;
			this->template store<true>((typename Q::value_t*)&q0, IConst<3>());
			q0.w = 0;
			return (q.inversion() * q0 * q).getVector();
		}
		this_t cross(const this_t& w) const {
			return base_t::I::Cross(base_t::m, w.m);
		}
		this_t verticalVector() const {
			auto ret = this_t(1,0,0) % *this;
			const value_t len_s = ret.len_sq();
			if(len_s < 1e-6f) {
				ret = this_t(0,0,1) % *this;
				return ret.normalization();
			}
			return ret * RSqrt(len_s);
		}
		//! 平面との交差点を算出
		template <bool A>
		auto planeDivide(const this_t& v, const PlaneT<value_t, A>& p) const {
			// 線分が平面をまたぐか
			const value_t distf = p.dot(*this);
			const value_t distb = p.dot(v);
			if(distf * distb >= 0)
				return std::make_tuple(this_t(), false);
		
			const value_t ratio = fabs(distf) / (fabs(distf) + fabs(distb));
			// 平面と線分の交点 -> tv
			this_t tv = v - (*this);
			tv *= ratio;
			const this_t cp = tv + (*this);
			return std::make_tuple(cp,true);
		}
		//! 平面にてベクトルを反転
		template <bool A>
		void flip(const PlaneT<value_t, A>& plane) {
			const auto& nml = plane.getNormal();
			const value_t d = plane.dot(*this);
			*this += nml*-d*2;
		}
	};
}
