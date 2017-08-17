//! 4次元ベクトル特有の関数など定義
//! (vector.hppからインクルード)
#pragma once

namespace frea {
	template <class R>
	struct wrap_spec<R,4> : wrap<R,4,wrap_spec<R,4>> {
		using base_t = wrap<R,4,wrap_spec<R,4>>;
		using base_t::base_t;
		wrap_spec() = default;
		wrap_spec(const base_t& b): base_t(b) {}

		auto asVec3Coord() const noexcept {
			auto tmp = *this;
			tmp.template makeEquality<3>();
			return (*this * base_t(base_t::I::Reciprocal(tmp))).template convert<3>();
		}
	};
	template <class T>
	struct tup_spec<T,4> : tup<T,4, tup_spec<T,4>> {
		using base_t = tup<T,4, tup_spec<T,4>>;
		using base_t::base_t;
		tup_spec() = default;
		tup_spec(const base_t& b): base_t(b) {}
		using this_t = tup_spec;
		using value_t = typename base_t::value_t;

		auto asVec3Coord() const noexcept {
			using v3 = typename base_t::template type_cn<3>;
			alignas(16) value_t tmp[4];
			this->template store<true>(tmp, lubee::IConst<3>());
			const value_t div = base_t::I::Reciprocal(tmp[3]);
			for(int i=0 ; i<3 ; i++)
				tmp[i] *= div;
			return v3(tmp, std::true_type());
		}
	};
}
