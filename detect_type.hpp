#pragma once
#include "type_fwd.hpp"
#include <type_traits>

namespace frea {
	template <class T, int N>
	std::true_type IsTuple(tup_spec<T,N>*);
	template <class T, int N, class S>
	std::true_type IsTuple(tup<T,N,S>*);
	std::false_type IsTuple(...);
	template <class T>
	using IsTuple_t = decltype(IsTuple((T*)nullptr));

	template <class T>
	struct is_vector : std::false_type {};
	template <class W, class D, int N>
	struct is_vector<VecT_spec<W,D,N>> : std::true_type {};
	template <class W, class D, class S>
	struct is_vector<VecT<W,D,S>> : std::true_type {};

	template <class T>
	struct is_wrap : std::false_type {};
	template <class R, int D, class S>
	struct is_wrap<wrap<R,D,S>> : std::true_type {};
	template <class R, int N>
	struct is_wrap<wrap_spec<R,N>> : std::true_type {};

	template <class T>
	struct is_matrix : std::false_type {};
	template <class V, int M, class S>
	struct is_matrix<MatT<V,M,S>> : std::true_type {};
	template <class V, int M, int N, class S>
	struct is_matrix<MatT_dspec<V,M,N,S>> : std::true_type {};
	template <class V, int M, int N>
	struct is_matrix<MatT_spec<V,M,N>> : std::true_type {};

	template <class T>
	struct is_wrapM : std::false_type {};
	template <class VW, int M, int N>
	struct is_wrapM<wrapM_spec<VW,M,N>> : std::true_type {};
	template <class VW, int M, class S>
	struct is_wrapM<wrapM<VW,M,S>> : std::true_type {};

	template <class T>
	struct is_quaternion : std::false_type {};
	template <class T, bool A>
	struct is_quaternion<QuatT<T,A>> : std::true_type {};

	template <class T>
	struct is_exp_quaternion : std::false_type {};
	template <class T, bool A>
	struct is_exp_quaternion<ExpQuatT<T,A>> : std::true_type {};

	template <class T>
	struct is_plane : std::false_type {};
	template <class T, bool A>
	struct is_plane<PlaneT<T,A>> : std::true_type {};
}
