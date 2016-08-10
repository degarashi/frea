#pragma once

namespace frea {
	#define countof(a) (sizeof(a)/sizeof(a[0]))
	#define DEF_OP(func, op) \
		template <class R> \
		static auto func(const R& t0, const R& t1) { \
			R ret; \
			for(int i=0 ; i<int(countof(t0.m)) ; i++) \
				ret.m[i] = t0.m[i] op t1.m[i]; \
			return ret; \
		}
	template <class T, bool F>
	struct _info {
		static T OneV() { return ~0; }
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
		static f_type OneV() {
			const i_type tmp(~0);
			return *(f_type*)&tmp;
		}
		#define DEF_OPF(func, op) \
			template <class R> \
			static auto func(const R& t0, const R& t1) { \
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

		static auto Zero() { return reg_t(0); }
		static auto One() { return reg_t(base_t::OneV()); }
		static auto LoadU(const value_t* v) { return reg_t(v); }
		static auto Load(const value_t* v) { return reg_t(v); }
		static auto Set1(const value_t& v) { return reg_t(v); }
		template <class... Ts>
		static auto Set(const Ts&... ts) { return reg_t(ts...); }
		template <int Pos>
		static auto MaskH(IConst<Pos>) { return reg_t(reg_t::TagMask, IConst<Pos>(), base_t::OneV()); }
		template <int Pos>
		static auto PickAt(IConst<Pos>) { return reg_t(reg_t::TagSet, IConst<Pos>(), base_t::OneV()); }
		template <int Pos>
		static value_t Pick(const reg_t& t) { return t.m[Pos]; }
		template <bool A2, int N2>
		static void Store(value_t* dst, const reg_t& t, BConst<A2>, IConst<N2>) {
			for(int i=0 ; i<N2+1 ; i++)
				*dst++ = t.m[i];
		}
		static value_t SumUp(const reg_t& t) {
			value_t sum = t.m[0];
			for(int i=1 ; i<capacity ; i++)
				sum += t.m[i];
			return sum;
		}
		template <class T2, bool A2>
		static reg_t Cast(const Data<T2,N,A2>& t) {
			reg_t ret;
			for(int i=0 ; i<capacity ; i++)
				ret.m[i] = *reinterpret_cast<const value_t*>(&t.m[i]);
			return ret;
		}
		template <class T2, bool A2>
		static reg_t Convert(const Data<T2,N,A2>& t) {
			reg_t ret;
			for(int i=0 ; i<capacity ; i++)
				ret.m[i] = t.m[i];
			return ret;
		}
	};
	#undef DEF_OP
}
