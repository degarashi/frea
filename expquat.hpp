#pragma once
#include "vector.hpp"
#include "angle.hpp"
#include "quaternion.hpp"
#include "fwd.hpp"

namespace frea {
	template <class T, bool A>
	struct ExpQuatT : Data<T,3,A>, lubee::op::Operator_Ne<ExpQuatT<T,A>> {
		using op_t = lubee::op::Operator_Ne<ExpQuatT<T,A>>;
		using base_t = Data<T,3,A>;
		using base_t::base_t;
		using rad_t = Radian<T>;
		using value_t = typename base_t::value_t;
		using vec_t = Vec_t<T,3,A>;
		using quat_t = QuatT<T,A>;

		ExpQuatT() = default;
		constexpr ExpQuatT(const quat_t& q) {
			constexpr auto Th = lubee::ThresholdF<value_t>(0.1);
			if(std::abs(q.w) >= 1.0-Th &&
				std::abs(q.x)+std::abs(q.y)+std::abs(q.z) <= Th)
			{
				// 無回転クォータニオンとみなす
				this->x = this->y = this->z = 0;
			} else {
				const value_t theta = std::acos(q.w);
				asVec3() = q.getAxis() * vec_t(theta);
			}
		}
		//! 異なる内部形式、アラインメントからの変換
		template <class T2, bool A2>
		constexpr ExpQuatT(const ExpQuatT<T2,A2>& q) noexcept:
			base_t(static_cast<const base_t&>(q))
		{}
		quat_t asQuat() const noexcept {
			const auto ret = getAngAxis();
			return quat_t::Rotation(ret.second, ret.first);
		}

		#define DEF_OP(op)	\
			ExpQuatT operator op (const ExpQuatT& e) const noexcept { \
				return asVec3() op e.asVec3(); \
			} \
			ExpQuatT operator op (const value_t& s) const noexcept { \
				return asVec3() op s; \
			} \
			using op_t::operator op;
		DEF_OP(+)
		DEF_OP(-)
		DEF_OP(*)
		DEF_OP(/)
		#undef DEF_OP

		bool operator == (const ExpQuatT& q) const noexcept {
			return asVec3() == q.asVec3();
		}
		value_t len_sq() const noexcept {
			return asVec3().len_sq();
		}
		value_t length() const noexcept {
			return asVec3().length();
		}
		vec_t& asVec3() noexcept {
			return reinterpret_cast<vec_t&>(*this);
		}
		const vec_t& asVec3() const noexcept {
			return reinterpret_cast<const vec_t&>(*this);
		}
		auto getAngAxis() const noexcept {
			std::pair<rad_t, vec_t> ret;
			ret.second = asVec3();
			const value_t theta = ret.second.length();		// = (angle/2)
			if(theta < lubee::ThresholdF<value_t>(0.1)) {
				// 無回転クォータニオンとする
				ret.first = rad_t(0);
				ret.second = vec_t(1,0,0);
			} else {
				ret.second /= theta;
				ret.first = rad_t(theta*2);
			}
			return ret;
		}
	};
	template <class T, bool A>
	inline std::ostream& operator << (std::ostream& os, const ExpQuatT<T,A>& q) {
		os << "ExpQuat: ";
		return q.asVec3().print(os);
	}
}
namespace std {
	template <class T, bool A>
	struct hash<frea::ExpQuatT<T,A>> {
		using eq_t = frea::ExpQuatT<T,A>;
		std::size_t operator()(const eq_t& q) const noexcept {
			using base_t = typename eq_t::base_t;
			return hash<base_t>()(static_cast<const base_t&>(q));
		}
	};
}
