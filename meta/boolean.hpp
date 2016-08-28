#pragma once
#include "constant_t.hpp"

namespace frea {
	template <class... Ts>
	struct And : std::true_type {};
	template <class T0, class... Ts>
	struct And<T0, Ts...> : BConst<T0{} && And<Ts...>{}> {};

	template <class... Ts>
	struct Or : std::true_type {};
	template <class T0, class... Ts>
	struct Or<T0, Ts...> : BConst<T0{} || Or<Ts...>{}> {};
}
