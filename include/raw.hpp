#pragma once

namespace frea {
	#define countof(a) (sizeof(a)/sizeof(a[0]))
	#define DEF_OP(func, op) \
		template <class R> \
		static auto func(const R& t0, const R& t1) noexcept { \
			R ret; \
			for(int i=0 ; i<int(countof(t0.m)) ; i++) \
				ret.m[i] = t0.m[i] op t1.m[i]; \
			return ret; \
		}
	template <class T, bool F>
	struct _info {
		static T OneV() noexcept { return ~0; }
		DEF_OP(And, &)
		DEF_OP(Or, |)
		DEF_OP(Xor, ^)
	};
	template <class T>
	struct _info<T, true> {
		static int32_t ToInt(float);
		static int64_t ToInt(double);
		template <class T2>
		using ToInt_t = decltype(ToInt(T2()));

		using i_type = _info<T,true>::ToInt_t<T>;
		using f_type = T;
		static f_type OneV() noexcept {
			const i_type tmp(~0);
			return *(f_type*)&tmp;
		}
		#define DEF_OPF(func, op) \
			template <class R> \
			static auto func(const R& t0, const R& t1) noexcept { \
				R ret; \
				for(int i=0 ; i<int(countof(t0.m)) ; i++) { \
					const i_type tmp = *reinterpret_cast<const i_type*>(&t0.m[i]) op \
										*reinterpret_cast<const i_type*>(&t1.m[i]); \
					ret.m[i] = *reinterpret_cast<const f_type*>(&tmp); \
				} \
				return ret; \
			}
		DEF_OPF(And, &)
		DEF_OPF(Or, |)
		DEF_OPF(Xor, ^)
		#undef DEF_OPF
	};
	template <class T, int N, bool A>
	struct info<Data<T,N,A>> : _info<T, std::is_floating_point<T>{}> {
		constexpr static int capacity = N;
		using base_t = _info<T, std::is_floating_point<T>{}>;
		using reg_t = Data<T,N,A>;
		using value_t = T;

		DEF_OP(Add, +)
		DEF_OP(Sub, -)
		DEF_OP(Mul, *)
		DEF_OP(Div, /)

		static auto Zero() noexcept { return reg_t(0); }
		static auto One() noexcept { return reg_t(base_t::OneV()); }
		static auto LoadU(const value_t* v) noexcept { return reg_t(v, std::false_type()); }
		static auto Load(const value_t* v) noexcept { return reg_t(v, std::true_type()); }
		static auto Set1(const value_t& v) noexcept { return reg_t(v); }
		template <std::size_t... Idx>
		static auto _Set(std::index_sequence<Idx...>, const value_t* ar) noexcept {
			return reg_t(ar[Idx]...);
		}
		template <class... Ts>
		static auto Set(const Ts&... ts) noexcept {
			// 一旦配列に格納して順序を逆転
			const value_t tmp[sizeof...(Ts)] = {static_cast<value_t>(ts)...};
			return _Set(seq::Reverse_t<std::make_index_sequence<sizeof...(Ts)>>(), tmp);
		}
		template <class... Ts>
		static auto SetR(const Ts&... ts) noexcept { return reg_t(ts...); }
		template <int Pos>
		static auto MaskH(IConst<Pos>) noexcept { return reg_t(TagMask, IConst<Pos>(), base_t::OneV()); }
		template <int Pos>
		static auto PickAt(IConst<Pos>) noexcept { return reg_t(TagSet, IConst<Pos>(), base_t::OneV()); }
		template <int Pos>
		static value_t Pick(const reg_t& t) noexcept { return t.m[Pos]; }
		template <bool A2, int N2>
		static void Store(value_t* dst, const reg_t& t, BConst<A2>, IConst<N2>) noexcept {
			static_assert(N2<N, "");
			for(int i=0 ; i<N2+1 ; i++)
				*dst++ = t.m[i];
		}
		static value_t SumUp(const reg_t& t) noexcept {
			value_t sum = t.m[0];
			for(int i=1 ; i<capacity ; i++)
				sum += t.m[i];
			return sum;
		}
		template <class T2, bool A2>
		static reg_t Cast(const Data<T2,N,A2>& t) noexcept {
			reg_t ret;
			for(int i=0 ; i<capacity ; i++)
				ret.m[i] = *reinterpret_cast<const value_t*>(&t.m[i]);
			return ret;
		}
		template <class T2, bool A2>
		static reg_t Convert(const Data<T2,N,A2>& t) noexcept {
			reg_t ret;
			for(int i=0 ; i<capacity ; i++)
				ret.m[i] = t.m[i];
			return ret;
		}
		static value_t GetMinValue(const reg_t& r) noexcept {
			value_t ret = r.m[0];
			for(int i=1 ; i<capacity ; i++)
				ret = std::min(ret, r.m[i]);
			return ret;
		}
		static value_t GetMaxValue(const reg_t& r) noexcept {
			value_t ret = r.m[0];
			for(int i=1 ; i<capacity ; i++)
				ret = std::max(ret, r.m[i]);
			return ret;
		}
		static bool Equal(const reg_t& r0, const reg_t& r1) noexcept {
			for(int i=0 ; i<capacity ; i++) {
				if(r0.m[i] != r1.m[i])
					return false;
			}
			return true;
		}
		static reg_t Min(const reg_t& m0, const reg_t& m1) noexcept {
			reg_t ret;
			for(int i=0 ; i<capacity ; i++)
				ret.m[i] = std::min(m0.m[i], m1.m[i]);
			return ret;
		}
		static reg_t Max(const reg_t& m0, const reg_t& m1) noexcept {
			reg_t ret;
			for(int i=0 ; i<capacity ; i++)
				ret.m[i] = std::max(m0.m[i], m1.m[i]);
			return ret;
		}
		static reg_t Absolute(const reg_t& r) noexcept {
			reg_t ret;
			for(int i=0 ; i<capacity ; i++)
				ret.m[i] = std::abs(r.m[i]);
			return ret;
		}
		static reg_t Lt(const reg_t& r0, const reg_t& r1) noexcept {
			reg_t ret;
			for(int i=0 ; i<capacity ; i++)
				ret.m[i] = (r0.m[i] < r1.m[i]) ? ~0 : 0;
			return ret;
		}
		static reg_t Sqrt(const reg_t& r) noexcept {
			reg_t ret;
			for(int i=0 ; i<capacity ; i++)
				ret.m[i] = std::sqrt(r.m[i]);
			return ret;
		}
		static bool IsNaN(const reg_t& r) noexcept {
			for(int i=0 ; i<capacity ; i++) {
				if(r.m[i] != r.m[i])
					return true;
			}
			return false;
		}
		static bool IsOutstanding(const reg_t& r) noexcept {
			for(int i=0 ; i<capacity ; i++) {
				const auto& v = r.m[i];
				if(v!=v || std::abs(v)>=std::numeric_limits<value_t>::infinity())
					return true;
			}
			return false;
		}
		static value_t Reciprocal(const value_t& v) noexcept {
			return 1 / v;
		}
		static reg_t Reciprocal(const reg_t& r) noexcept {
			reg_t ret;
			for(int i=0 ; i<capacity ; i++)
				ret.m[i] = 1/r.m[i];
			return ret;
		}
		static reg_t Cross(const reg_t& a, const reg_t& b) noexcept {
			reg_t c;
			c.m[0] = a.m[1]*b.m[2] - a.m[2]*b.m[1];
			c.m[1] = a.m[2]*b.m[0] - a.m[0]*b.m[2];
			c.m[2] = a.m[0]*b.m[1] - a.m[1]*b.m[0];
			return c;
		}
	};
	#undef DEF_OP
}
