#pragma once
#include <utility>

namespace frea {
	template <class T, class... Args>
	class IsCallable {
		private:
			template <class F>
			static std::true_type _check(decltype(std::declval<F>()(std::declval<Args>()...), (void)0)*);
			template <class F>
			static std::false_type _check(...);
		public:
			using type = decltype(_check<T>(nullptr));
			constexpr static bool value = type::value;
	};
	template <class T, class... Args>
	using IsCallable_t = typename IsCallable<T, Args...>::type;
}
