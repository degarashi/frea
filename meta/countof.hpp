#pragma once
#include <cstdint>

namespace frea {
	template <std::size_t N>
	struct GetCountOf_helper {
		using type = char [N];
	};

	template <class T, std::size_t N>
	typename GetCountOf_helper<N>::type& GetCountOf(T (&)[N]);
	//! 配列の要素数を取得する (配列でない場合はエラー)
	#define countof(e)		sizeof(::spn::GetCountOf(e))

	template <class T>
	char GetCountOfNA(T);
	template <class T, std::size_t N>
	typename GetCountOf_helper<N>::type& GetCountOfNA(T (&)[N]);
	//! 配列の要素数を取得する (配列でない場合は1が返る)
	#define countof_na(e)	sizeof(::frea::GetCountOfNA(e))
}
