//! 2次元ベクトル特有の関数など定義
//! (vector.hppからインクルード)
#pragma once

namespace frea {
	template <class R>
	struct wrap_spec<R,2> : wrap<R,2, wrap_spec<R,2>> {
		using base_t = wrap<R, 2, wrap_spec<R,2>>;
		using base_t::base_t;
		using value_t = typename base_t::value_t;
		using this_t = wrap_spec;
		wrap_spec() = default;
		wrap_spec(const base_t& b): base_t(b) {}

		struct D2 {
			value_t x,y;
		};
		D2 asD2() const noexcept {
			D2 ret;
			this->template store<false>(reinterpret_cast<value_t*>(&ret), lubee::IConst<1>());
			return ret;
		}
		//! counter clockwise
		/*! 半時計回りが正数を返す */
		value_t ccw(const this_t& v) const noexcept {
			const auto a = asD2(),
					 b = v.asD2();
			return a.x*b.y - a.y*b.x;
		}
		//! clockwise
		/*! 時計回りが正数を返す */
		value_t cw(const this_t& v) const noexcept {
			const auto a = asD2(),
					b = v.asD2();
			return -a.x*b.y + a.y*b.x;
		}
	};
}
