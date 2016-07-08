#pragma once
#include "constant_t.hpp"

namespace frea {
	//! クラスがwidthフィールドを持っていればintegral_constantでそれを返す
	template <class T>
	IConst<T::width> _GetWidthT(decltype(T::width)*);
	template <class T>
	IConst<0> _GetWidthT(...);
	template <class T>
	decltype(_GetWidthT<T>(nullptr)) GetWidthT();
	//! クラスがheightフィールドを持っていればintegral_constantでそれを返す
	template <class T>
	IConst<T::height> _GetHeightT(decltype(T::height)*);
	template <class T>
	IConst<0> _GetHeightT(...);
	template <class T>
	decltype(_GetHeightT<T>(nullptr)) GetHeightT();
	//! クラスがwidthとheightフィールドを持っていればintegral_pairでそれを返す
	template <class T>
	IPair<T::height, T::width> _GetWidthHeightT(decltype(T::width)*, decltype(T::height)*);
	template <class T>
	IPair<0, T::width> _GetWidthHeightT(decltype(T::width)*, ...);
	template <class T>
	IPair<0,0> _GetWidthHeightT(...);
	template <class T>
	decltype(_GetWidthHeightT<T>(nullptr,nullptr)) GetWidthHeightT();
}
