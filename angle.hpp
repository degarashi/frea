#pragma once
#include "constant.hpp"
#include "lubee/meta/enable_if.hpp"
#include "lubee/meta/check_macro.hpp"
#include "lubee/range.hpp"
#include "lubee/tostring.hpp"
#include <cmath>

namespace frea {
	// Degree角度を示すタグ構造体
	struct Degree_t {};
	// Radian角度を示すタグ構造体
	struct Radian_t {};

	//! 角度変換ヘルパークラス
	template <class From, class To>
	struct ConvertAngle;
	// 同じ単位なので変換の必要なし
	template <class TAG>
	struct ConvertAngle<TAG, TAG> {
		template <class T>
		T operator()(const T ang) const noexcept {
			return ang;
		}
	};
	template <>
	struct ConvertAngle<Degree_t, Radian_t> {
		template <class T>
		T operator()(const T deg) const noexcept {
			return deg / 180 * Pi<T>;
		}
	};
	template <>
	struct ConvertAngle<Radian_t, Degree_t> {
		template <class T>
		T operator()(const T rad) const noexcept {
			return rad / Pi<T> * 180;
		}
	};

	template <class TAG>
	struct AngleInfo;
	template <>
	struct AngleInfo<Degree_t> {
		template <class T>
		constexpr static T one_rotation{360};
		template <class Dummy=void>
		constexpr static char name[]= {"Degree"},
							name_short[] = {"Deg"};
	};
	template <class Dummy> const char AngleInfo<Degree_t>::name[];
	template <class Dummy> const char AngleInfo<Degree_t>::name_short[];
	template <>
	struct AngleInfo<Radian_t> {
		template <class T>
		constexpr static T one_rotation{2*Pi<T>};
		template <class Dummy=void>
		constexpr static char name[] = {"Radian"},
							name_short[] = {"Rad"};
	};
	template <class Dummy> const char AngleInfo<Radian_t>::name[];
	template <class Dummy> const char AngleInfo<Radian_t>::name_short[];

	//! 角度クラス
	/*! 混同を避けるために数値への変換は明示的
		他形式との変換は暗黙的 */
	template <class TAG, class V>
	class Angle {
		public:
			using value_t = V;
		private:
			using tag_type = TAG;
			value_t	_angle;
			template <class V2>
			static auto _Mod(const V2& f, const V2& m, std::false_type) noexcept {
				return std::fmod(f, m);
			}
			template <class V2>
			static auto _Mod(const V2& f, const V2& m, std::true_type) noexcept {
				return f % m;
			}
			template <class V2>
			static auto Mod(const V2& f, const V2& m) noexcept {
				return _Mod(f, m, lubee::HasOp_Mod_t<V2,V2>());
			}

		public:
			constexpr static value_t OneRotationAng = AngleInfo<tag_type>::template one_rotation<value_t>;
			constexpr static lubee::Range<value_t> OneRotationRange = {0, OneRotationAng},
											HalfRotationRange = {-OneRotationAng/2, OneRotationAng/2};

			Angle() = default;
			//! 他タイプの角度で角度指定
			template <class TAG2, class V2>
			constexpr Angle(const Angle<TAG2,V2>& ang) noexcept {
				_angle = ConvertAngle<TAG2, TAG>()(ang.get());
			}
			//! floatで直接角度指定
			explicit constexpr Angle(const value_t& ang) noexcept:
				_angle(ang)
			{}
			constexpr static Angle Rotation(value_t r) noexcept {
				return Angle(OneRotationAng * r);
			}
			Angle& operator = (const Angle& ang) = default;
			template <class TAG2, class V2>
			Angle& operator = (const Angle<TAG2,V2>& ang) noexcept {
				_angle = Angle(ang).get();
				return *this;
			}

			// 浮動少数点数型なら明示的に角度を出力
			template <class T2, ENABLE_IF((std::is_floating_point<T2>{}))>
			explicit operator T2() const noexcept {
				return _angle;
			}
			template <class TAG2, class V2=value_t>
			constexpr Angle<TAG2,V2> convert() const noexcept {
				return Angle<TAG2,V2>(*this);
			}
			void set(const value_t ang) noexcept {
				_angle = ang;
			}
			constexpr value_t get() const noexcept {
				return _angle;
			}
			Angle distance(const Angle& r) const noexcept {
				return Angle(std::abs(get() - r.get()));
			}
			void semicircle() noexcept {
				constexpr auto OR = OneRotationAng;
				auto ang = Mod(_angle, OR);
				if(ang > OR/2)
					ang -= OR;
				_angle = ang;
			}
			//! 角度をループ(0から2Piの間に収める)
			void single() noexcept {
				constexpr auto OR = OneRotationAng;
				auto ang = Mod(_angle, OR);
				if(ang < 0)
					ang += OneRotationAng;
				_angle = ang;
			}
			//! 角度を一定の範囲に制限(直の値)
			void rangeValue(const lubee::Range<value_t>& r) noexcept {
				if(_angle < r.from)
					_angle = r.from;
				else if(_angle > r.to)
					_angle = r.to;
			}
			//! 角度を一定の範囲に制限(直の値)
			void range(const lubee::Range<Angle>& r) noexcept {
				rangeValue({r.from.get(), r.to.get()});
			}

			template <class TAG2, class V2>
			Angle operator + (const Angle<TAG2,V2>& ang) const noexcept {
				return Angle(_angle + Angle(ang).get());
			}
			template <class TAG2, class V2>
			Angle& operator += (const Angle<TAG2,V2>& ang) noexcept {
				_angle += Angle(ang).get();
				return *this;
			}
			template <class TAG2, class V2>
			Angle operator - (const Angle<TAG2,V2>& ang) const noexcept {
				return Angle(_angle - Angle(ang).get());
			}
			template <class TAG2, class V2>
			Angle& operator -= (const Angle<TAG2,V2>& ang) noexcept {
				_angle -= Angle(ang).get();
				return *this;
			}
			Angle operator * (const value_t r) const noexcept {
				return Angle(_angle * r);
			}
			Angle& operator *= (const value_t r) noexcept {
				_angle *= r;
				return *this;
			}
			Angle operator / (const value_t r) const noexcept {
				return Angle(_angle / r);
			}
			Angle& operator /= (const value_t r) noexcept {
				_angle /= r;
				return *this;
			}
			Angle operator -() const noexcept {
				return Angle(-_angle);
			}

			// ---- 比較演算子の定義 ----
			#define DEF_OP(op) \
				template <class TAG2, class V2> \
				bool operator op (const Angle<TAG2,V2>& ang) const noexcept { \
					return _angle op Angle(ang).get(); \
				}
			DEF_OP(==)
			DEF_OP(!=)
			DEF_OP(<)
			DEF_OP(<=)
			DEF_OP(>)
			DEF_OP(>=)
			#undef DEF_OP
			template <class Ar>
			value_t save_minimal(const Ar&) const {
				return _angle;
			}
			template <class Ar>
			void load_minimal(const Ar&, const value_t& v) {
				_angle = v;
			}

			// -------- Luaへのエクスポート用 --------
			using Deg_t = Angle<Degree_t, value_t>;
			using Rad_t = Angle<Radian_t, value_t>;
			Angle luaAddD(const Deg_t& a) const {
				return *this + a;
			}
			Angle luaAddR(const Rad_t& a) const {
				return *this + a;
			}
			Angle luaSubD(const Deg_t& a) const {
				return *this - a;
			}
			Angle luaSubR(const Rad_t& a) const {
				return *this - a;
			}
			Deg_t luaToDegree() const {
				return *this;
			}
			Rad_t luaToRadian() const {
				return *this;
			}
			Angle luaMulF(const float s) const {
				return *this * s;
			}
			Angle luaDivF(const float s) const {
				return *this / s;
			}
			Angle luaInvert() const {
				return -*this;
			}
			bool luaLessthan(const Angle& a) const {
				return *this < a;
			}
			bool luaLessequal(const Angle& a) const {
				return *this <= a;
			}
			bool luaEqual(const Angle& a) const {
				return *this == a;
			}
			std::string luaToString() const {
				return lubee::ToString(*this);
			}
	};
	template <class TAG, class V>
	constexpr V Angle<TAG,V>::OneRotationAng;
	template <class TAG, class V>
	constexpr lubee::Range<V> Angle<TAG,V>::OneRotationRange;
	template <class TAG, class V>
	constexpr lubee::Range<V> Angle<TAG,V>::HalfRotationRange;

	template <class T>
	using Degree = Angle<Degree_t, T>;
	using DegF = Degree<float>;
	using DegD = Degree<double>;
	template <class T>
	using Radian = Angle<Radian_t, T>;
	using RadF = Radian<float>;
	using RadD = Radian<double>;

	template <class TAG, class V>
	inline std::ostream& operator << (std::ostream& os, const Angle<TAG,V>& ang) {
		return os << ang.get() << "(" << AngleInfo<TAG>::template name_short<> << ")";
	}
}
namespace std {
	template <class TAG, class V>
	struct hash<frea::Angle<TAG,V>> {
		std::size_t operator()(const frea::Angle<TAG,V>& r) const {
			const auto a = r.get();
			return std::hash<std::decay_t<decltype(a)>>()(a);
		}
	};
}
