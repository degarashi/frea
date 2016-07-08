#pragma once
#include "vector.hpp"
#include "angle.hpp"
#include "quaternion.hpp"

namespace frea {
	template <class T, bool A>
	struct ExpQuatT : Data<T,3,A> {
		using base_t = Data<T,3,A>;
		using base_t::base_t;
		using rad_t = Radian<T>;
		using value_t = typename base_t::value_t;
		using vec_t = Vec_t<T,3,A>;
		using quat_t = QuatT<T,A>;

		ExpQuatT() = default;
		constexpr ExpQuatT(const quat_t& q) noexcept {
			if(std::fabs(q.w) >= 1.0-1e-6) {
				// 無回転クォータニオンとみなす
				this->x = this->y = this->z = 0;
			} else {
				const value_t theta = std::acos(q.w);
				asVec3() = q.getAxis() * vec_t(theta);
			}
		}
		template <class T2, bool A2>
		constexpr ExpQuatT(const ExpQuatT<T2,A2>& q) noexcept:
			base_t(static_cast<const base_t&>(q))
		{}
		quat_t asQuat() const {
			const auto ret = getAngAxis();
			return quat_t::Rotation(ret.second, ret.first);
		}

		#define DEF_OP(op)	\
			ExpQuatT&& operator op (const ExpQuatT& q) && { \
				return std::move(*this op##= q); \
			} \
			ExpQuatT operator op (const ExpQuatT& q) const& { \
				return asVec3() * q.asVec3(); \
			} \
			ExpQuatT&& operator op (const value_t& s) && { \
				return std::move(*this op##= s); \
			} \
			ExpQuatT operator op (const value_t& s) const& { \
				return asVec3() * s; \
			} \
			template <class T2> \
			ExpQuatT& operator op##= (const T2& t) { \
				return *this = *this op t; \
			} \
			template <class T2> \
			ExpQuatT&& operator op##= (const T2& t) && { \
				return std::move(*this = *this op t); \
			}
		DEF_OP(+)
		DEF_OP(-)
		DEF_OP(*)
		DEF_OP(/)
		#undef DEF_OP

		value_t len_sq() const {
			return asVec3().len_sq();
		}
		value_t length() const {
			return asVec3().length();
		}
		vec_t& asVec3() noexcept {
			return reinterpret_cast<vec_t&>(*this);
		}
		const vec_t& asVec3() const noexcept {
			return static_cast<const vec_t&>(*this);
		}
		std::pair<rad_t,vec_t> getAngAxis() const {
			auto axis = asVec3();
			const value_t theta = axis.length();		// = (angle/2)
			if(theta < 1e-6) {
				// 無回転クォータニオンとする
				return std::make_pair(rad_t(0), vec_t(1,0,0));
			}
			axis /= theta;
			return std::make_pair(rad_t(theta*2), axis);
		}
	};
}
