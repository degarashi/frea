#pragma once

namespace frea {
	template <class T, int N, class S>
	struct tup;
	template <class T, int N>
	struct tup_spec;

	template <class R, int D, class S>
	struct wrap;
	template <class R, int N>
	struct wrap_spec;
	template <class W, class D, class S>
	struct VecT;
	template <class W, class D, int N>
	struct VecT_spec;

	template <class VW, int M, class S>
	class wrapM;
	template <class VW, int M, int N>
	struct wrapM_spec;
	template <class V, int M, class S>
	struct MatT;
	template <class V, int M, int N, class S>
	struct MatT_dspec;
	template <class V, int M, int N>
	struct MatT_spec;

	template <class T, bool A>
	struct QuatT;
	template <class T, bool A>
	struct ExpQuatT;
	template <class T, bool A>
	struct PlaneT;
}
