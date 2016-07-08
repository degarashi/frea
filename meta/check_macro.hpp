#pragma once
#include <type_traits>

namespace frea {
	template <class T, class Idx>
	std::false_type HasIndex(...);
	template <class T, class Idx>
	std::true_type HasIndex(std::decay_t<decltype(std::declval<T>()[std::declval<Idx>()])>*);
	template <class T, class Idx>
	using HasIndex_t = decltype(HasIndex<T,Idx>(nullptr));

	#define DEF_HasOp(name, op) \
		template <class T0, class T1> \
		std::false_type HasOp_##name(...); \
		template <class T0, class T1> \
		std::true_type HasOp_##name(std::decay_t<decltype(std::declval<T0>() op std::declval<T1>())>*); \
		template <class T0, class T1> \
		using HasOp_##name##_t = decltype(HasOp_##name<T0,T1>(nullptr));
	DEF_HasOp(Plus, +)
	DEF_HasOp(Minus, -)
	DEF_HasOp(Mul, *)
	DEF_HasOp(Mod, %)
	#undef DEF_HasOp
}

//! 特定のメソッドを持っているかチェック
/*! DEF_HASMETHOD(method)
	HasMethod_method<class_name>(nullptr) -> std::true_type or std::false_type */
#define DEF_HASMETHOD_BASE(method) \
	template <class T> \
	std::false_type HasMethod_##method(...); \
	template <class T> \
	using HasMethod_##method##_t = decltype(HasMethod_##method<T>(nullptr));
#define DEF_HASMETHOD(method) \
	template <class T> \
	std::true_type HasMethod_##method(decltype(&T::method) = &T::method); \
	DEF_HASMETHOD_BASE(method)
#define DEF_HASMETHOD_T(method) \
	template <class T, class... Ts> \
	std::true_type HasMethod_##method(decltype(&T::template method<Ts...>)); \
	DEF_HASMETHOD_BASE(method)

//! 特定のオーバーロードされたメソッドを持っているかチェック
/*! DEF_HASMETHOD_OV(name, method)
	予めクラスにfrea::none method(...); を持たせておく
	暗黙の変換も含まれる
	HasMethod_name<class_name>() -> std::true_type or std::false_type */
#define DEF_HASMETHOD_OV(name, ...) \
	template <class T, class... Args> \
	auto CheckMethod_##name() -> std::integral_constant<bool, !std::is_same<::frea::none_t, decltype(std::declval<T>().method(std::declval<Args>()...))>::value>; \
	template <class T> \
	auto HasMethod_##name() -> decltype(CheckMethod_##name<T, __VA_ARGS__>()); \
	template <class T> \
	using HasMethod_##name##_t = decltype(HasMethod_##method<T>(nullptr));

//! クラスが特定の名前の定数値を持っているかチェック
/*! DEF_HASTYPE(type_name)
	HasType_type_name<class_name>(nullptr) -> std::true_type or std::false_type */
#define DEF_HASTYPE(name) \
	template <class T> \
	std::true_type HasType_##name(decltype(T::name)*); \
	template <class T> \
	std::false_type HasType_##name(...); \
	template <class T> \
	using HasType_##name##_t = decltype(HasType_##name<T>(nullptr));

//! クラスが特定の名前のtype aliasを持っているかチェック
/*! DEF_HASTYPE_T(type_name)
	HasTypeT_type_name<class_name>(nullptr) -> std::true_type or std::false_type */
#define DEF_HASTYPE_T(name) \
	template <class T> \
	std::true_type HasTypeT_##name(typename T::type*); \
	template <class T> \
	std::false_type HasTypeT_##name(...); \
	template <class T> \
	using HasTypeT_##name##_t = decltype(HasTypeT_##name<T>(nullptr));


