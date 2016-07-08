#pragma once
#include <type_traits>

namespace frea {
	//! 任意のバイト数でアラインメント
	template <std::size_t N>
	struct alignas(N) Align {};
	template <bool A>
	struct Align16 : std::conditional_t<A, Align<16>, Align<4>> {};
}
