#pragma once
#include "constant_t.hpp"
#include <functional>

namespace frea {
	//! std::tupleのハッシュを計算(適当実装)
	struct TupleHash {
		template <class Tup>
		static std::size_t get(std::size_t value, const Tup& /*tup*/, IConst<-1>) { return value; }
		template <class Tup, int N>
		static std::size_t get(std::size_t value, const Tup& tup, IConst<N>) {
			const auto& t = std::get<N>(tup);
			std::size_t h = std::hash<std::decay_t<decltype(t)>>()(t);
			value = (value ^ (h<<(h&0x07))) ^ (h>>3);
			return get(value, tup, IConst<N-1>());
		}
		template <class... Ts>
		std::size_t operator()(const std::tuple<Ts...>& tup) const {
			return get(0xdeadbeef * 0xdeadbeef, tup, IConst<sizeof...(Ts)-1>());
		}
	};
}
