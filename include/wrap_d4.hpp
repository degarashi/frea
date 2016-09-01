//! 4次元ベクトル特有の関数など定義
//! (vector.hppからインクルード)
#pragma once

namespace frea {
	template <class R>
	struct wrap_spec<R,4> : wrap<R,4,wrap_spec<R,4>> {
		using base_t = wrap<R,4,wrap_spec<R,4>>;
		using base_t::base_t;

		typename base_t::template type_cn<3> asVec3Coord() const {
			auto tmp = *this;
			tmp.template makeEquality<3>();
			return this->template convert<3>() * base_t::I::Reciprocal(tmp);
		}
	};
	template <class T>
	struct tup_spec<T,4> : tup<T,4, tup_spec<T,4>> {
		using base_t = tup<T,4, tup_spec<T,4>>;
		using base_t::base_t;
		using this_t = tup_spec;
		using value_t = typename base_t::value_t;

		auto asVec3Coord() const {
			using v3 = typename base_t::template type_cn<3>;
			alignas(16) value_t tmp[4];
			this->template store<true>(tmp, IConst<3>());
			const value_t div = base_t::I::Reciprocal(tmp[3]);
			for(int i=0 ; i<3 ; i++)
				tmp[i] *= div;
			return v3(tmp, std::true_type());
		}
	};
}
