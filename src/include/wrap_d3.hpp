//! 3次元ベクトル特有の関数など定義
//! (vector.hppからインクルード)
#pragma once
#include <tuple>

namespace frea {
	template <class W>
	W VerticalVector(const W& w) noexcept {
		auto ret = W(1,0,0) % w;
		const auto len_s = ret.len_sq();
		if(len_s < 1e-6) {
			ret = W(0,0,1) % w;
			return ret.normalization();
		}
		return ret / std::sqrt(len_s);
	}
	template <class R>
	struct wrap_spec<R,3> : wrap<R,3,wrap_spec<R,3>> {
		using base_t = wrap<R,3,wrap_spec<R,3>>;
		using base_t::base_t;
		wrap_spec() = default;
		wrap_spec(const base_t& b): base_t(b) {}
		using this_t = wrap_spec;
		using value_t = typename base_t::value_t;

		using base_t::operator *;
		template <class Q,
				 ENABLE_IF(is_quaternion<Q>{})>
		auto operator * (const Q& q) const noexcept ->
			std::decay_t<decltype(std::declval<QuatT<typename Q::value_t, true>>().getVector())>
		{
			QuatT<typename Q::value_t, true> q0;
			this->template store<true>((typename Q::value_t*)&q0, lubee::IConst<2>());
			q0.w = 0;
			return (q * q0 * q.inversion()).getVector();
		}
		this_t cross(const this_t& w) const noexcept {
			return base_t::I::Cross(base_t::m, w.m);
		}
		this_t operator % (const this_t& w) const noexcept {
			return cross(w);
		}
		this_t verticalVector() const noexcept {
			return VerticalVector(*this);
		}
	};
	template <class T>
	struct tup_spec<T,3> : tup<T,3, tup_spec<T,3>> {
		using base_t = tup<T,3, tup_spec<T,3>>;
		using base_t::base_t;
		tup_spec() = default;
		tup_spec(const base_t& b): base_t(b) {}
		using this_t = tup_spec;

		using base_t::operator *;
		template <class Q,
				 ENABLE_IF(is_quaternion<Q>{})>
		auto operator * (const Q& q) const noexcept ->
			std::decay_t<decltype(std::declval<QuatT<typename Q::value_t, true>>().getVector())>
		{
			QuatT<typename Q::value_t, true> q0;
			this->template store<true>(reinterpret_cast<typename Q::value_t*>(&q0), lubee::IConst<2>());
			q0.w = 0;
			return (q * q0 * q.inversion()).getVector();
		}
		using value_t = typename base_t::value_t;
		this_t cross(const this_t& t) const noexcept {
			alignas(16) value_t a[3],
								b[3],
								c[3];
			this->template store<true>(a, lubee::IConst<2>());
			t.template store<true>(b, lubee::IConst<2>());
			c[0] = a[1]*b[2] - a[2]*b[1];
			c[1] = a[2]*b[0] - a[0]*b[2];
			c[2] = a[0]*b[1] - a[1]*b[0];
			return this_t(c, lubee::BConst<true>());
		}
		this_t operator % (const this_t& w) const noexcept {
			return cross(w);
		}
		this_t verticalVector() const noexcept {
			return VerticalVector(*this);
		}
	};
}
