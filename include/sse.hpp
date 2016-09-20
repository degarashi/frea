#pragma once
#include <xmmintrin.h>
#include <emmintrin.h>

#if defined(SSE)
#if SSE <= 2
	#define SUMVEC(r)	{ auto tmp = _mm_shuffle_ps(r, r, _MM_SHUFFLE(0,1,2,3)); \
		tmp = Add(r, tmp); \
		r = _mm_shuffle_ps(tmp,tmp, _MM_SHUFFLE(1,0,0,1)); \
		r = Add(r, tmp); }
#else
	// rの各要素を足し合わせる
	#define SUMVEC(r)	{ r = _mm_hadd_ps(r,r); \
		r = _mm_hadd_ps(r,r); }
#endif
namespace frea {
	template <>
	struct info<__m128> {
		constexpr static int capacity = 4;
		using reg_t = __m128;
		using regI_t = __m128i;
		using regD_t = __m128d;
		using value_t = float;
		constexpr static auto Add = &_mm_add_ps,
								Sub = &_mm_sub_ps,
								Mul = &_mm_mul_ps,
								Div = &_mm_div_ps,
								And = &_mm_and_ps,
								Or = &_mm_or_ps,
								Xor = &_mm_xor_ps,
								Min = &_mm_min_ps,
								Max = &_mm_max_ps,
								Lt = &_mm_cmplt_ps;
		constexpr static auto Sqrt = &_mm_sqrt_ps;
		constexpr static auto Set1 = &_mm_set1_ps;
		constexpr static auto Set = &_mm_set_ps,
							SetR = &_mm_setr_ps;
		constexpr static auto Zero = &_mm_setzero_ps;
		constexpr static auto Load = &_mm_load_ps,
							LoadU = &_mm_loadu_ps;

		#define AsReg(w,z,y,x)	_mm_castsi128_ps(_mm_set_epi32(w,z,y,x))
		static auto AbsMask() noexcept {
			constexpr uint32_t m = 0x7fffffff;
			return AsReg(m,m,m,m);
		}
		static auto One() noexcept { return AsReg(-1,-1,-1,-1); }
		static auto MaskH(lubee::IConst<0>) noexcept { return AsReg(0,0,0,-1); }
		static auto MaskH(lubee::IConst<1>) noexcept { return AsReg(0,0,-1,-1); }
		static auto MaskH(lubee::IConst<2>) noexcept { return AsReg(0,-1,-1,-1); }
		static auto MaskH(lubee::IConst<3>) noexcept { return AsReg(-1,-1,-1,-1); }
		static auto PickAt(lubee::IConst<0>) noexcept { return AsReg(0,0,0,-1); }
		static auto PickAt(lubee::IConst<1>) noexcept { return AsReg(0,0,-1,0); }
		static auto PickAt(lubee::IConst<2>) noexcept { return AsReg(0,-1,0,0); }
		static auto PickAt(lubee::IConst<3>) noexcept { return AsReg(-1,0,0,0); }
		#undef AsReg

		template <int N>
		static value_t Pick(const reg_t& t) noexcept {
			static_assert(N>=0 && N<capacity, "invalid index");
			const reg_t t2 = _mm_shuffle_ps(t, t, _MM_SHUFFLE(N,N,N,N));
			value_t ret;
			_mm_store_ss(&ret, t2);
			return ret;
		}
		template <bool A>
		static void Store(value_t* dst, const reg_t& t, lubee::BConst<A>, lubee::IConst<0>) noexcept {
			_mm_store_ss(dst, t);
		}
		template <bool A>
		static void Store(value_t* dst, const reg_t& t, lubee::BConst<A>, lubee::IConst<1>) noexcept {
			_mm_storel_pi(reinterpret_cast<__m64*>(dst), t);
		}
		template <bool A>
		static void Store(value_t* dst, const reg_t& t, lubee::BConst<A>, lubee::IConst<2>) noexcept {
			_mm_storel_pi(reinterpret_cast<__m64*>(dst), t);
			_mm_store_ss(dst+2, _mm_shuffle_ps(t, t, _MM_SHUFFLE(2,2,2,2)));
		}
		static void Store(value_t* dst, const reg_t& t, std::false_type, lubee::IConst<3>) noexcept {
			_mm_storeu_ps(dst, t);
		}
		static void Store(value_t* dst, const reg_t& t, std::true_type, lubee::IConst<3>) noexcept {
			_mm_store_ps(dst, t);
		}
		static value_t SumUp(const reg_t& r) noexcept {
			auto tr = r;
			SUMVEC(tr)
			value_t tmp;
			Store(&tmp, tr, std::false_type(), lubee::IConst<0>());
			return tmp;
		}
		static reg_t Cast(const regI_t& t) noexcept {
			return _mm_castsi128_ps(t);
		}
		static reg_t Cast(const regD_t& t) noexcept {
			return _mm_castpd_ps(t);
		}
		static reg_t Convert(const regI_t& t) noexcept {
			return _mm_cvtepi32_ps(t);
		}
		static reg_t Convert(const regD_t& t) noexcept {
			return _mm_cvtpd_ps(t);
		}
		static void Transpose(reg_t& r0, reg_t& r1, reg_t& r2, reg_t& r3) noexcept {
			reg_t xmt[4];
			xmt[0] = _mm_unpacklo_ps(r0, r2);
			xmt[1] = _mm_unpacklo_ps(r1, r3);
			xmt[2] = _mm_unpackhi_ps(r0, r2);
			xmt[3] = _mm_unpackhi_ps(r1, r3);
			r0 = _mm_unpacklo_ps(xmt[0], xmt[1]);
			r1 = _mm_unpackhi_ps(xmt[0], xmt[1]);
			r2 = _mm_unpacklo_ps(xmt[2], xmt[3]);
			r3 = _mm_unpackhi_ps(xmt[2], xmt[3]);
		}
		static void Transpose(reg_t& r0, reg_t& r1, reg_t& r2) noexcept {
			reg_t r3 = Zero();
			Transpose(r0, r1, r2, r3);
		}
		static bool Equal(const reg_t& r0, const reg_t& r1) noexcept {
			auto t0 = _mm_cmpeq_ps(r0, r1);
			t0 = _mm_and_ps(t0, _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(1,0,3,2)));
			t0 = _mm_and_ps(t0, _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(0,1,2,3)));
			return _mm_cvttss_si32(t0) != 0;
		}
		static reg_t _NaNChk(const reg_t& r0) noexcept {
			const auto res = Or(_mm_cmple_ps(r0, Zero()),
							_mm_cmpgt_ps(r0, Zero()));
			return _mm_andnot_ps(res, One());
		}
		static bool IsNaN(const reg_t& r0) noexcept {
			const auto res = _NaNChk(r0);
			return SumUp(res) != 0;
		}
		static bool IsOutstanding(const reg_t& r) noexcept {
			auto r0 = And(r, AbsMask());
			const auto r_inf = And(Set1(std::numeric_limits<value_t>::infinity()), AbsMask());
			r0 = _mm_cmpeq_ps(r0, r_inf);
			r0 = Or(r0, _NaNChk(r));
			return SumUp(r0) != 0;
		}
		static value_t Reciprocal(const value_t& v) noexcept {
			const reg_t r = Reciprocal(Set1(v));
			value_t ret;
			Store(&ret, r, lubee::BConst<false>(), lubee::IConst<0>());
			return ret;
		}
		static reg_t Reciprocal(reg_t r) noexcept {
			reg_t tmp(_mm_rcp_ps(r));
			r = Mul(r, Mul(tmp,tmp));
			tmp = Add(tmp, tmp);
			return Sub(tmp, r);
		}
		static reg_t Cross(const reg_t& t0, const reg_t& t1) noexcept {
			auto r0 = t0,
				 r1 = t1;
			// r0[y,z,x]
			auto m0 = _mm_shuffle_ps(r0, r0, _MM_SHUFFLE(0,0,2,1)),
			// r0[z,x,y]
					m1 = _mm_shuffle_ps(r0, r0, _MM_SHUFFLE(0,1,0,2)),
			// r1[z,x,y]
					m2 = _mm_shuffle_ps(r1, r1, _MM_SHUFFLE(0,1,0,2)),
			// r1[y,z,x]
					m3 = _mm_shuffle_ps(r1, r1, _MM_SHUFFLE(0,0,2,1));
			r0 = Mul(m0,m2);
			r1 = Mul(m1,m3);
			r0 = Sub(r0, r1);
			return r0;
		}
		static reg_t Absolute(const reg_t& r) noexcept {
			return And(AbsMask(), r);
		}
		template <class F>
		static value_t _GetValue(const F func, const reg_t& r) noexcept {
			auto tmp = _mm_shuffle_ps(r, r, _MM_SHUFFLE(0,0,2,3));
			tmp = func(r, tmp);
			auto tmp2 = _mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(0,0,0,1));
			tmp = func(tmp, tmp2);
			value_t ret;
			Store(&ret, tmp, std::false_type(), lubee::IConst<0>());
			return ret;
		}
		static value_t GetMinValue(const reg_t& r) noexcept {
			return _GetValue(Min, r);
		}
		static value_t GetMaxValue(const reg_t& r) noexcept {
			return _GetValue(Max, r);
		}
	};
	template <>
	struct info<__m128i> {
		constexpr static int capacity = 4;
		using reg_t = __m128i;
		using regF_t = __m128;
		using regD_t = __m128d;
		using value_t = int32_t;

		static reg_t Mul(reg_t m0, reg_t m1) noexcept {
			reg_t	s0 = _mm_shuffle_epi32(m0, _MM_SHUFFLE(0,1,0,0)),
					s1 = _mm_shuffle_epi32(m1, _MM_SHUFFLE(0,1,0,0));
			reg_t t0 = _mm_mul_epu32(s0, s1);
			t0 = _mm_shuffle_epi32(t0, _MM_SHUFFLE(0,0,2,0));
			s0 = _mm_shuffle_epi32(m0, _MM_SHUFFLE(0,3,0,2));
			s1 = _mm_shuffle_epi32(m1, _MM_SHUFFLE(0,3,0,2));
			reg_t t1 = _mm_mul_epu32(s0, s1);
			t1 = _mm_shuffle_epi32(t1, _MM_SHUFFLE(0,0,2,0));
			return _mm_movelh_ps(t0, t1);
		}
		static reg_t Div(reg_t m0, reg_t m1) noexcept {
			__m128 t0 = _mm_cvtepi32_ps(m0),
					t1 = _mm_cvtepi32_ps(m1);
			t0 = _mm_div_ps(t0, t1);
			return _mm_cvttps_epi32(t0);
		}
		static reg_t Load(const value_t* i) noexcept {
			return _mm_set_epi32(i[3],i[2],i[1],i[0]);
		}
		static reg_t Min(reg_t m0, reg_t m1) noexcept {
			return _Proc<&_mm_cmplt_epi32>(m0, m1);
		}
		static reg_t Max(reg_t m0, reg_t m1) noexcept {
			return _Proc<&_mm_cmpgt_epi32>(m0, m1);
		}
		constexpr static auto Add = &_mm_add_epi32,
								Sub = &_mm_sub_epi32,
								And = &_mm_and_si128,
								Or = &_mm_or_si128,
								Xor = &_mm_xor_si128,
								Lt = &_mm_cmplt_epi32;
		constexpr static auto LoadU = &Load;
		constexpr static auto Set1 = &_mm_set1_epi32;
		constexpr static auto Set = &_mm_set_epi32,
							SetR = &_mm_setr_epi32;
		constexpr static auto Zero = &_mm_setzero_si128;

		#define AsReg(w,z,y,x)	_mm_set_epi32(w,z,y,x)
		static auto One() noexcept { return AsReg(-1,-1,-1,-1); }
		static auto MaskH(lubee::IConst<0>) noexcept { return AsReg(0,0,0,-1); }
		static auto MaskH(lubee::IConst<1>) noexcept { return AsReg(0,0,-1,-1); }
		static auto MaskH(lubee::IConst<2>) noexcept { return AsReg(0,-1,-1,-1); }
		static auto MaskH(lubee::IConst<3>) noexcept { return AsReg(-1,-1,-1,-1); }
		static auto PickAt(lubee::IConst<0>) noexcept { return AsReg(0,0,0,-1); }
		static auto PickAt(lubee::IConst<1>) noexcept { return AsReg(0,0,-1,0); }
		static auto PickAt(lubee::IConst<2>) noexcept { return AsReg(0,-1,0,0); }
		static auto PickAt(lubee::IConst<3>) noexcept { return AsReg(-1,0,0,0); }
		#undef AsReg

		template <__m128i (*Cmp)(reg_t, reg_t)>
		static reg_t _Proc(reg_t m0, reg_t m1) noexcept {
			auto mask = Cmp(m0, m1);
			m0 = And(m0, mask);
			mask = Xor(mask, One());
			m1 = And(m1, mask);
			return Or(m0, m1);
		}

		template <bool A>
		static void Store(value_t* dst, const reg_t& t, lubee::BConst<A>, lubee::IConst<0>) noexcept {
			alignas(16) value_t tmp[capacity];
			_mm_store_si128(reinterpret_cast<reg_t*>(tmp), t);
			dst[0] = tmp[0];
		}
		template <bool A>
		static void Store(value_t* dst, const reg_t& t, lubee::BConst<A>, lubee::IConst<1>) noexcept {
			alignas(16) value_t tmp[capacity];
			_mm_store_si128(reinterpret_cast<reg_t*>(tmp), t);
			dst[0] = tmp[0];
			dst[1] = tmp[1];
		}
		template <bool A>
		static void Store(value_t* dst, const reg_t& t, lubee::BConst<A>, lubee::IConst<2>) noexcept {
			alignas(16) value_t tmp[capacity];
			_mm_store_si128(reinterpret_cast<reg_t*>(tmp), t);
			dst[0] = tmp[0];
			dst[1] = tmp[1];
			dst[2] = tmp[2];
		}
		static void Store(value_t* dst, const reg_t& t, std::false_type, lubee::IConst<3>) noexcept {
			_mm_storeu_si128(reinterpret_cast<reg_t*>(dst), t);
		}
		static void Store(value_t* dst, const reg_t& t, std::true_type, lubee::IConst<3>) noexcept {
			_mm_store_si128(reinterpret_cast<reg_t*>(dst), t);
		}
		static value_t SumUp(const reg_t& r) noexcept {
			auto tr = r;
			SUMVEC(tr)
			value_t tmp;
			Store(&tmp, tr, std::false_type(), lubee::IConst<0>());
			return tmp;
		}
		static reg_t Cast(const regF_t& t) noexcept {
			return _mm_castps_si128(t);
		}
		static reg_t Cast(const regD_t& t) noexcept {
			return _mm_castpd_si128(t);
		}
		static reg_t Convert(const regF_t& t) noexcept {
			return _mm_cvtps_epi32(t);
		}
		static reg_t Convert(const regD_t& t) noexcept {
			return _mm_cvtpd_epi32(t);
		}
		static bool Equal(const reg_t& r0, const reg_t& r1) noexcept {
			auto t0 = _mm_cmpeq_epi32(r0, r1);
			t0 = _mm_and_ps(t0, _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(1,0,3,2)));
			t0 = _mm_and_ps(t0, _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(0,1,2,3)));
			return _mm_cvttss_si32(t0) != 0;
		}
		static void Transpose(reg_t& r0, reg_t& r1, reg_t& r2, reg_t& r3) noexcept {
			using F = info<__m128>;
			F::Transpose(
				reinterpret_cast<regF_t&>(r0),
				reinterpret_cast<regF_t&>(r1),
				reinterpret_cast<regF_t&>(r2),
				reinterpret_cast<regF_t&>(r3)
			);
		}
		static void Transpose(reg_t& r0, reg_t& r1, reg_t& r2) noexcept {
			using F = info<__m128>;
			F::Transpose(
				reinterpret_cast<regF_t&>(r0),
				reinterpret_cast<regF_t&>(r1),
				reinterpret_cast<regF_t&>(r2)
			);
		}
		template <int N>
		static value_t Pick(const reg_t& t) noexcept {
			static_assert(N>=0 && N<capacity, "invalid index");
			const reg_t t2 = _mm_shuffle_epi32(t, _MM_SHUFFLE(N,N,N,N));
			value_t ret;
			Store(&ret, t2, lubee::BConst<false>(), lubee::IConst<1>());
			return ret;
		}
		static bool IsNaN(const reg_t&) noexcept {
			return false;
		}
		static bool IsOutstanding(const reg_t&) noexcept {
			return false;
		}
		static reg_t Absolute(const reg_t& r) noexcept {
			const auto lt = Lt(r, Zero());
			auto tmp = Sub(r, And(lt, _mm_set_epi32(1,1,1,1)));
			return Xor(tmp, And(lt, One()));
		}
		template <class F>
		static value_t _GetValue(const F func, const reg_t& r) noexcept {
			auto tmp = _mm_shuffle_epi32(r, _MM_SHUFFLE(0,0,2,3));
			tmp = func(r, tmp);
			auto tmp2 = _mm_shuffle_epi32(tmp, _MM_SHUFFLE(0,0,0,1));
			tmp = func(tmp, tmp2);
			value_t ret;
			Store(&ret, tmp, std::false_type(), lubee::IConst<0>());
			return ret;
		}
		static value_t GetMinValue(const reg_t& r) noexcept {
			return _GetValue(Min, r);
		}
		static value_t GetMaxValue(const reg_t& r) noexcept {
			return _GetValue(Max, r);
		}
	};
	template <>
	struct info<__m128d> {
		constexpr static int capacity = 2;
		using reg_t = __m128d;
		using regF_t = __m128;
		using regI_t = __m128i;
		using value_t = double;
		constexpr static auto Add = &_mm_add_pd,
								Sub = &_mm_sub_pd,
								Mul = &_mm_mul_pd,
								Div = &_mm_div_pd,
								And = &_mm_and_pd,
								Or = &_mm_or_pd,
								Xor = &_mm_xor_pd,
								Min = &_mm_min_pd,
								Max = &_mm_max_pd,
								Lt = &_mm_cmplt_pd;
		constexpr static auto Set1 = &_mm_set1_pd;
		constexpr static auto Set = &_mm_set_pd,
							SetR = &_mm_setr_pd;
		constexpr static auto Zero = &_mm_setzero_pd;
		constexpr static auto Load = &_mm_load_pd,
							LoadU = &_mm_loadu_pd;

		#define AsReg(w,z,y,x)	_mm_castsi128_pd(_mm_set_epi32(w,z,y,x))
		static auto AbsMask() noexcept {
			constexpr uint32_t m7f = 0x7fffffff,
								mff = 0xffffffff;
			return AsReg(m7f, mff, m7f, mff);
		}
		static auto One() noexcept { return AsReg(-1,-1,-1,-1); }
		static auto MaskH(lubee::IConst<0>) noexcept { return AsReg(0,0,-1,-1); }
		static auto MaskH(lubee::IConst<1>) noexcept { return AsReg(-1,-1,-1,-1); }
		static auto PickAt(lubee::IConst<0>) noexcept { return AsReg(0,0,-1,-1); }
		static auto PickAt(lubee::IConst<1>) noexcept { return AsReg(-1,-1,0,0); }
		#undef AsReg
		template <bool A>
		static void Store(value_t* dst, const reg_t& t, lubee::BConst<A>, lubee::IConst<0>) noexcept {
			_mm_store_sd(dst, t);
		}
		static void Store(value_t* dst, const reg_t& t, std::true_type, lubee::IConst<1>) noexcept {
			_mm_store_pd(dst, t);
		}
		static void Store(value_t* dst, const reg_t& t, std::false_type, lubee::IConst<1>) noexcept {
			_mm_storeu_pd(dst, t);
		}
		static value_t SumUp(const reg_t& r) noexcept {
			auto tr = _mm_shuffle_pd(r, r, 0b11);
			tr = Add(tr, r);
			value_t ret;
			Store(&ret, tr, std::false_type(), lubee::IConst<0>());
			return ret;
		}

		static reg_t Cast(const regI_t& t) noexcept {
			return _mm_castsi128_pd(t);
		}
		static reg_t Cast(const regF_t& t) noexcept {
			return _mm_castps_pd(t);
		}
		static reg_t Convert(const regI_t& t) noexcept {
			return _mm_cvtepi32_pd(t);
		}
		static reg_t Convert(const regF_t& t) noexcept {
			return _mm_cvtps_pd(t);
		}
		static bool Equal(const reg_t& r0, const reg_t& r1) noexcept {
			auto t0 = _mm_cmpeq_pd(r0, r1);
			t0 = _mm_and_ps(t0, _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(1,0,3,2)));
			t0 = _mm_and_ps(t0, _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(0,1,2,3)));
			return _mm_cvttss_si32(t0) != 0;
		}
		template <int N>
		static value_t Pick(const reg_t& t) noexcept {
			static_assert(N>=0 && N<capacity, "invalid index");
			const reg_t t2 = _mm_shuffle_pd(t, t, N);
			value_t ret;
			Store(&ret, t2, lubee::BConst<false>(), lubee::IConst<0>());
			return ret;
		}
		static bool IsNaN(const reg_t& r0) noexcept {
			auto r_zero = Zero();
			auto res = Or(_mm_cmple_pd(r0, r_zero),
							_mm_cmpgt_pd(r0, r_zero));
			res = _mm_andnot_pd(res, One());
			const value_t f = SumUp(res);
			return f != 0;
		}
		static bool IsOutstanding(const reg_t& r) noexcept {
			const auto f = std::numeric_limits<value_t>::infinity();
			auto r_inf = _mm_load1_pd(&f);
			auto r0 = And(r, AbsMask());
			auto r1 = Or(_mm_cmple_pd(r0, Zero()),
						_mm_cmpgt_pd(r0, Zero()));
			r0 = _mm_cmpeq_pd(r0, r_inf);

			r1 = _mm_andnot_pd(r1, One());
			r0 = Or(r0, r1);
			return SumUp(r0) != 0;
		}
		static value_t Reciprocal(const value_t& v) noexcept {
			return 1.0 / v;
		}
		static reg_t Reciprocal(reg_t r) noexcept {
			return Div(Set1(1), r);
		}
		static reg_t Absolute(const reg_t& r) noexcept {
			return And(AbsMask(), r);
		}
		template <class F>
		static value_t _GetValue(const F func, const reg_t& r) noexcept {
			auto tmp = _mm_shuffle_pd(r, r, 0b11);
			tmp = func(r, tmp);
			value_t ret;
			Store(&ret, tmp, std::false_type(), lubee::IConst<0>());
			return ret;
		}
		static value_t GetMinValue(const reg_t& r) noexcept {
			return _GetValue(Min, r);
		}
		static value_t GetMaxValue(const reg_t& r) noexcept {
			return _GetValue(Max, r);
		}
	};
	template <int N, bool A> SVec_t<__m128, N, A> info_detect(float);
	template <int N, bool A> SVec_t<__m128, N, A> info_detect(__m128);

	template <int N, bool A> SVec_t<__m128d, N, A> info_detect(double);
	template <int N, bool A> SVec_t<__m128d, N, A> info_detect(__m128d);

	template <int N, bool A> SVec_t<__m128i, N, A> info_detect(int32_t);
	template <int N, bool A> SVec_t<__m128i, N, A> info_detect(__m128i);
}
#endif
