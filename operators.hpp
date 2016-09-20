#pragma once
#include <utility>
#include "lubee/meta/enable_if.hpp"

namespace frea {
	namespace op {
		#define _NOEXCEPT(e) noexcept((Self&)std::declval<Self>() e (T&)std::declval<T>())
		#define NOEXCEPT(e) noexcept(_NOEXCEPT(e))
		#define NOEXCEPT_2(e, e2) noexcept(_NOEXCEPT(e) && _NOEXCEPT(e2))
		#define DEF_OP(op) \
			template <class T> \
			auto& operator op##= (const T& t) & NOEXCEPT(op) { \
				auto& self = static_cast<Self&>(*this); \
				return self = self op t; \
			} \
			template <class T> \
			auto&& operator op##= (const T& t) && NOEXCEPT(op##=) { \
				return std::move(*this op##= t); \
			} \
			template <class T> \
			decltype(auto) operator op (const T& t) && { \
				return static_cast<Self&>(*this) op t; \
			}
			// template <class T, ENABLE_IF((std::is_assignable<Self,decltype((Self&)std::declval<Self>() op (T&)std::declval<T>())>::value))> \
			// decltype(auto) operator op (const T& t) && { \
			// 	return std::move(static_cast<Self&>(*this) op##= t); \
			// } \
			// template <class T, ENABLE_IF(!(std::is_assignable<Self,decltype((Self&)std::declval<Self>() op (T&)std::declval<T>())>::value))> \
			// decltype(auto) operator op (const T& t) && { \
			// 	return static_cast<Self&>(*this) op t; \
			// }
		template <class Self>
		struct PlusMinus {
			DEF_OP(+)
			DEF_OP(-)
		};
		template <class Self>
		struct Mul {
			DEF_OP(*)
		};
		template <class Self>
		struct Div {
			DEF_OP(/)
		};
		template <class Self>
		struct MulDiv : Mul<Self>, Div<Self> {};
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
		struct Ne {
			template <class T>
			bool operator != (const T& t) const NOEXCEPT(==) {
				const auto& self = static_cast<const Self&>(*this);
				return !(self == t);
			}
		};
		template <class Self>
		struct Ge {
			template <class T>
			bool operator > (const T& t) const NOEXCEPT_2(<, ==) {
				const auto& self = static_cast<const Self&>(*this);
				return !(self < t || self == t);
			}
			template <class T>
			bool operator >= (const T& t) const NOEXCEPT(<) {
				const auto& self = static_cast<const Self&>(*this);
				return !(self < t);
			}
			template <class T>
			bool operator <= (const T& t) const NOEXCEPT(<) {
				const auto& self = static_cast<const Self&>(*this);
				return (self < t) || (self == t);
			}
		};
		template <class Self>
		struct Cmp : Ne<Self>, Ge<Self> {};
		#undef _NOEXCEPT
		#undef NOEXCEPT
		#undef NOEXCEPT_2
		template <class Self>
		struct Operator : Arithmetic<Self>, Logical<Self>, Cmp<Self> {};
		template <class Self>
		struct Operator_Ne : Arithmetic<Self>, Logical<Self>, Ne<Self> {};
	}
}
