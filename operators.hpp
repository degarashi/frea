#pragma once
#include <utility>

namespace frea {
	namespace op {
		#define NOEXCEPT(e) //noexcept(noexcept((Self&)std::declval<Self>() e (T&)std::declval<T>()))
		#define NOEXCEPT_R(e) //noexcept(noexcept(std::declval<Self>() e std::declval<T>()))
		#define DEF_OP(op) \
			template <class T> \
			auto& operator op##= (const T& t) & NOEXCEPT(op##=) { \
				auto& self = static_cast<Self&>(*this); \
				return self = self op t; \
			} \
			template <class T> \
			auto&& operator op##= (const T& t) && NOEXCEPT_R(op##=) { \
				return std::move(*this op##= t); \
			} \
			template <class T> \
			auto&& operator op (const T& t) && NOEXCEPT_R(op) { \
				return std::move(static_cast<Self&>(*this) op##= t); \
			}
		template <class Self>
		struct PlusMinus {
			DEF_OP(+)
			DEF_OP(-)
		};
		template <class Self>
		struct MulDiv {
			DEF_OP(*)
			DEF_OP(/)
		};
		template <class Self>
		struct Arithmetic : PlusMinus<Self>, MulDiv<Self> {};
		template <class Self>
		struct Logical {
			DEF_OP(&)
			DEF_OP(|)
			DEF_OP(^)
		};
		#undef DEF_OP
		template <class Self>
		struct Compare {
			template <class T>
			bool operator != (const T& t) const NOEXCEPT(!=) {
				const auto& self = static_cast<const Self&>(*this);
				return !(self == t);
			}
			template <class T>
			bool operator > (const T& t) const NOEXCEPT(>) {
				const auto& self = static_cast<const Self&>(*this);
				return !(self < t || self == t);
			}
			template <class T>
			bool operator >= (const T& t) const NOEXCEPT(>=) {
				const auto& self = static_cast<const Self&>(*this);
				return !(self < t);
			}
			template <class T>
			bool operator <= (const T& t) const NOEXCEPT(<=) {
				const auto& self = static_cast<const Self&>(*this);
				return (self < t) || (self == t);
			}
		};
		#undef NOEXCEPT
		template <class Self>
		struct Operator : Arithmetic<Self>, Logical<Self>, Compare<Self> {};
	}
}
