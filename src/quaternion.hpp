#pragma once
#include "matrix.hpp"
#include "lubee/src/error.hpp"
#include "compare.hpp"
#include "lubee/src/ieee754.hpp"
#include "exception.hpp"
#include "lubee/src/compare.hpp"
#include "fwd.hpp"

namespace frea {
	template <class T, bool A>
	struct QuatT : Data<T,4,A>, lubee::op::Operator_Ne<QuatT<T,A>> {
		using op_t = lubee::op::Operator_Ne<QuatT<T,A>>;
		using base_t = Data<T,4,A>;
		using base_t::base_t;
		using vec_t = Vec_t<T,3,A>;
		using vec4_t = Vec_t<T,4,A>;
		using rad_t = Radian<T>;
		using value_t = typename base_t::value_t;
		using mat3_t = Mat_t<T,3,3, A>;
		using mat4_t = Mat_t<T,4,4, A>;
		using exp_t = ExpQuatT<value_t, A>;

		constexpr static T ZeroLen_Th = lubee::ThresholdF<value_t>(0.2),
							Theta_Th = lubee::ThresholdF<value_t>(0.56);
		QuatT() = default;
		operator const base_t&() = delete;
		// 違う要素Quatからの変換
		template <class T2, bool A2>
		constexpr QuatT(const QuatT<T2,A2>& q) noexcept:
			base_t(static_cast<const base_t&>(q))
		{}
		constexpr static QuatT Identity() noexcept {
			return {0,0,0,1};
		}
		constexpr QuatT(const vec_t& v, const value_t& w) noexcept:
			base_t{v.x, v.y, v.z, w}
		{}
		static QuatT FromAxisF(const value_t& f11, const value_t& f12, const value_t& f13,
								const value_t& f21, const value_t& f22, const value_t& f23,
								const value_t& f31, const value_t& f32, const value_t& f33) noexcept
		{
			// 最大成分を検索
			const value_t elem[4] = {f11 - f22 - f33 + 1,
									-f11 + f22 - f33 + 1,
									-f11 - f22 + f33 + 1,
									f11 + f22 + f33 + 1};
			int idx = 0;
			for(int i=1 ; i<4 ; i++) {
				if(elem[i] > elem[idx])
					idx = i;
			}
			D_Expect(elem[idx] >= 0, "invalid matrix error")

			// 最大要素の値を算出
			QuatT res;
			auto *const qr = res.m;
			const value_t v = std::sqrt(elem[idx]) * 0.5;
			qr[idx] = v;
			value_t mult = 0.25 / v;
			switch(idx) {
				case 0: // x
					qr[1] = (f12 + f21) * mult;
					qr[2] = (f31 + f13) * mult;
					qr[3] = (f23 - f32) * mult;
					break;
				case 1: // y
					qr[0] = (f12 + f21) * mult;
					qr[2] = (f23 + f32) * mult;
					qr[3] = (f31 - f13) * mult;
					break;
				case 2: // z
					qr[0] = (f31 + f13) * mult;
					qr[1] = (f23 + f32) * mult;
					qr[3] = (f12 - f21) * mult;
					break;
				case 3: // w
					qr[0] = (f23 - f32) * mult;
					qr[1] = (f31 - f13) * mult;
					qr[2] = (f12 - f21) * mult;
					break;
			}
			return res;
		}
		static QuatT FromMat(const mat3_t& m) noexcept {
			return FromAxis(m.m[0], m.m[1], m.m[2]);
		}
		static QuatT FromAxis(const vec_t& xA,
								const vec_t& yA,
								const vec_t& zA) noexcept
		{
			return FromAxisF(xA.x, xA.y, xA.z,
							yA.x, yA.y, yA.z,
							zA.x, zA.y, zA.z);
		}
		static QuatT FromMatAxis(const vec_t& xA,
								const vec_t& yA,
								const vec_t& zA) noexcept
		{
			return FromAxisF(xA.x, yA.x, zA.x,
							xA.y, yA.y, zA.y,
							xA.z, yA.z, zA.z);
		}
		static QuatT RotationYPR(const rad_t yaw, const rad_t pitch, const rad_t roll) noexcept {
			QuatT q = Identity();
			// roll
			q.rotateZ(roll);
			// pitch
			q.rotateX(-pitch);
			// yaw
			q.rotateY(yaw);
			return q;
		}
		static QuatT RotationX(rad_t ang) noexcept {
			ang *= 0.5;
			return QuatT(std::sin(ang.get()), 0, 0, std::cos(ang.get()));
		}
		static QuatT RotationY(rad_t ang) noexcept {
			ang *= 0.5;
			return QuatT(0, std::sin(ang.get()), 0, std::cos(ang.get()));
		}
		static QuatT RotationZ(rad_t ang) noexcept {
			ang *= 0.5;
			return QuatT(0, 0, std::sin(ang.get()), std::cos(ang.get()));
		}
		static QuatT Rotation(const vec_t& axis, const rad_t ang) noexcept {
			D_Expect(std::abs(1-axis.length()) < lubee::ThresholdF<value_t>(0.4), "invalid axis")
			D_Expect(std::abs(ang.get()) < rad_t::OneRotationAng, "invalid angle range")
			const auto C = std::cos(ang.get()/2),
						S = std::sin(ang.get()/2);
			const vec_t taxis = axis * S;
			return QuatT(taxis, C);
		}
		static QuatT LookAt(const vec_t& dir, const vec_t& up) noexcept {
			vec_t t_up = up;
			vec_t rv = t_up % dir;
			const value_t len_s = rv.len_sq();
			if(len_s < ZeroLen_Th) {
				// 真上か真下を向いている
				// upベクトルは適当に定める
				t_up = dir.verticalVector();
				rv = t_up % dir;
			} else {
				rv /= std::sqrt(len_s);
				t_up = dir % rv;
			}
			return FromAxis(rv, t_up, dir);
		}
		static QuatT Rotation(const vec_t& from, const vec_t& to) noexcept {
			vec_t rAxis = from % to;
			if(rAxis.len_sq() < ZeroLen_Th)
				return QuatT::Identity();
			rAxis.normalize();
			const value_t d = std::acos(lubee::Saturate<value_t>(from.dot(to), 1.0));
			return Rotation(rAxis, rad_t(d));
		}
		static QuatT SetLookAt(const Axis::e targetAxis, const Axis::e baseAxis, const vec_t& baseVec, const vec_t& at, const vec_t& pos) {
			if(targetAxis == baseAxis)
				throw InvalidAxis();
			// [0] = target
			// [1] = right
			// [2] = base
			int axF[3] = {targetAxis, 0, baseAxis};
			// axF[1]へtarget, base以外の軸フラグを設定
			switch((1<<targetAxis) | (1<<baseAxis)) {
				// 110
				case 0x06:
					axF[1] = Axis::X;
					break;
				// 101
				case 0x05:
					axF[1] = Axis::Y;
					break;
				// 011
				case 0x03:
					axF[1] = Axis::Z;
					break;
				default:
					D_Expect(false, "invalid axis flag")
			}
			vec_t axis[3];
			auto &vTarget = axis[axF[0]],
				&vOther = axis[axF[1]],
				&vBase = axis[axF[2]];
			vTarget = at - pos;
			if(vTarget.normalize() < ZeroLen_Th)
				throw NoValidAxis();
			vOther = baseVec.cross(vTarget);
			if(vOther.normalize() < ZeroLen_Th)
				throw NoValidAxis();
			vBase = vTarget.cross(vOther);
			const value_t b = (int(targetAxis) == (int(baseAxis)+1)%3) ? 1 : -1;
			vOther *= b;
			return FromAxis(axis[0], axis[1], axis[2]);
		}
		void rotateX(const rad_t ang) noexcept {
			*this = RotationX(ang) * *this;
		}
		void rotateY(const rad_t ang) noexcept {
			*this = RotationY(ang) * *this;
		}
		void rotateZ(const rad_t ang) noexcept {
			*this = RotationZ(ang) * *this;
		}
		void rotate(const vec_t& axis, const rad_t ang) noexcept {
			*this = Rotation(axis, ang) * *this;
		}
		void identity() noexcept {
			*this = Identity();
		}
		void normalize() noexcept {
			*this = normalization();
		}
		void conjugate() noexcept {
			*this = conjugation();
		}
		void invert() noexcept {
			*this = inversion();
		}
		QuatT scale(const value_t& s) const noexcept {
			return Identity().slerp(*this, s);
		}

		#define DEF_SCALAR(op) \
			QuatT operator op (const value_t& v) const noexcept { \
				return asVec4() op v; \
			} \
			using op_t::operator op;
		DEF_SCALAR(+)
		DEF_SCALAR(-)
		DEF_SCALAR(*)
		DEF_SCALAR(/)
		#undef DEF_SCALAR

		#define DEF_QUAT(op) \
			QuatT operator op (const QuatT& v) const noexcept { \
				return asVec4() op v.asVec4(); \
			}
		DEF_QUAT(+)
		DEF_QUAT(-)
		#undef DEF_QUAT

		QuatT operator * (const QuatT& q) const noexcept {
			return {
				this->w*q.x + this->x*q.w + this->y*q.z - this->z*q.y,
				this->w*q.y - this->x*q.z + this->y*q.w + this->z*q.x,
				this->w*q.z + this->x*q.y - this->y*q.x + this->z*q.w,
				this->w*q.w - this->x*q.x - this->y*q.y - this->z*q.z
			};
		}

		QuatT rotationX(const rad_t ang) const noexcept {
			return RotationX(ang) * *this;
		}
		QuatT rotationY(const rad_t ang) const noexcept {
			return RotationY(ang) * *this;
		}
		QuatT rotationZ(const rad_t ang) const noexcept {
			return RotationZ(ang) * *this;
		}
		QuatT rotation(const vec_t& axis, rad_t ang) const noexcept {
			return Rotation(axis, ang) * *this;
		}
		const vec4_t& asVec4() const noexcept {
			return reinterpret_cast<const vec4_t&>(*this);
		}
		QuatT normalization() const noexcept {
			return asVec4().normalization();
		}
		QuatT conjugation() const noexcept {
			return
				QuatT(
					-this->x,
					-this->y,
					-this->z,
					this->w
				);
		}
		QuatT inversion() const noexcept {
			return conjugation() / len_sq();
		}
		value_t len_sq() const noexcept {
			return asVec4().len_sq();
		}
		value_t length() const noexcept {
			return asVec4().length();
		}
		rad_t angle() const noexcept {
			return rad_t(std::acos(lubee::Saturate<value_t>(this->w, 1.0))*2);
		}
		const vec_t& getVector() const noexcept {
			return reinterpret_cast<const vec_t&>(*this);
		}
		vec_t getAxis() const {
			auto s_theta = std::sqrt(1.0 - Square(this->w));
			if(s_theta < ZeroLen_Th)
				throw NoValidAxis();
			s_theta = 1.0 / s_theta;
			return vec_t(this->x*s_theta, this->y*s_theta, this->z*s_theta);
		}

		value_t dot(const QuatT& q) const noexcept {
			return asVec4().dot(q.asVec4());
		}
		bool operator == (const QuatT& q) const noexcept {
			return asVec4() == q.asVec4();
		}
		//! 球面線形補間
		QuatT slerp(const QuatT& q, const value_t& t) const noexcept {
			const auto ac = lubee::Saturate<value_t>(dot(q), 0.0, 1.0);
			const auto theta = std::acos(ac),
						S = std::sin(theta);
			if(std::abs(S) < Theta_Th)
				return *this;
			QuatT rq = *this * (std::sin(theta*(1-t)) / S);
			rq += q * (std::sin(theta * t) / S);
			return rq;
		}

		#define ELEM00 (1-2*this->y*this->y-2*this->z*this->z)
		#define ELEM01 (2*this->x*this->y+2*this->w*this->z)
		#define ELEM02 (2*this->x*this->z-2*this->w*this->y)
		#define ELEM10 (2*this->x*this->y-2*this->w*this->z)
		#define ELEM11 (1-2*this->x*this->x-2*this->z*this->z)
		#define ELEM12 (2*this->y*this->z+2*this->w*this->x)
		#define ELEM20 (2*this->x*this->z+2*this->w*this->y)
		#define ELEM21 (2*this->y*this->z-2*this->w*this->x)
		#define ELEM22 (1-2*this->x*this->x-2*this->y*this->y)
		//! 回転を行列表現した時のX軸
		vec_t getXAxis() const noexcept {
			return {ELEM00, ELEM10, ELEM20};
		}
		//! 正規直行座標に回転を掛けた後のX軸
		vec_t getXAxisInv() const noexcept {
			return {ELEM00, ELEM01, ELEM02};
		}
		//!< 回転を行列表現した時のY軸
		vec_t getYAxis() const noexcept {
			return {ELEM01, ELEM11, ELEM21};
		}
		//!< 正規直行座標に回転を掛けた後のY軸
		vec_t getYAxisInv() const noexcept {
			return {ELEM10, ELEM11, ELEM12};
		}
		//!< 回転を行列表現した時のZ軸
		vec_t getZAxis() const noexcept {
			return {ELEM02, ELEM12, ELEM22};
		}
		//!< 正規直行座標に回転を掛けた後のZ軸
		vec_t getZAxisInv() const noexcept {
			return {ELEM20, ELEM21, ELEM22};
		}
		//!< X軸に回転を適用したベクトル
		vec_t getRight() const noexcept {
			return {ELEM00, ELEM01, ELEM02};
		}
		//!< Y軸に回転を適用したベクトル
		vec_t getUp() const noexcept {
			return {ELEM10, ELEM11, ELEM12};
		}
		//!< Z軸に回転を適用したベクトル
		vec_t getDir() const noexcept {
			return {ELEM20, ELEM21, ELEM22};
		}
		//! 行列変換(3x3)
		mat3_t asMat33() const noexcept {
			return {
				ELEM00, ELEM01, ELEM02,
				ELEM10, ELEM11, ELEM12,
				ELEM20, ELEM21, ELEM22
			};
		}
		//! 行列変換(4x4)
		mat4_t asMat44() const noexcept {
			return {
				ELEM00, ELEM01, ELEM02, 0,
				ELEM10, ELEM11, ELEM12, 0,
				ELEM20, ELEM21, ELEM22, 0,
				0,0,0,1
			};
		}
		#undef ELEM00
		#undef ELEM01
		#undef ELEM02
		#undef ELEM10
		#undef ELEM11
		#undef ELEM12
		#undef ELEM20
		#undef ELEM21
		#undef ELEM22
		//! (デバッグ用)2つのクォータニオンにどれだけ差があるかの目安
		value_t distance(const QuatT& q) const noexcept {
			// 単純に値の差分をとる
			const value_t f[4] = {
				this->x - q.x,
				this->y - q.y,
				this->z - q.z,
				this->w - q.w
			};
			value_t sum = 0;
			for(int i=0 ; i<4 ; i++)
				sum += std::abs(f[i]);
			return sum;
		}
		exp_t asExpQuat() const {
			return exp_t(*this);
		}
		// -------- Luaへのエクスポート用 --------
		static QuatT Lua_Rotation(const Vec3& axis, const RadF ang) {
			return Rotation(axis, ang);
		}
		static QuatT Lua_RotationFromTo(const Vec3& from, const Vec3& to) {
			return Rotation(from, to);
		}
		QuatT luaRotation(const vec_t& axis, const rad_t ang) const noexcept {
			return Rotation(axis, ang);
		}
		QuatT luaRotationFromTo(const vec_t& from, const vec_t& to) const noexcept {
			return Rotation(from, to);
		}
		QuatT luaAddQ(const QuatT& q) const noexcept {
			return *this + q;
		}
		QuatT luaSubQ(const QuatT& q) const noexcept {
			return *this - q;
		}
		QuatT luaMulQ(const QuatT& q) const noexcept {
			return *this * q;
		}
		QuatT luaMulF(const float s) const noexcept {
			return *this * s;
		}
		QuatT luaDivF(const float s) const noexcept {
			return *this / s;
		}
		bool luaEqual(const QuatT& q) const noexcept {
			return *this == q;
		}
		std::string luaToString() const {
			return lubee::ToString(*this);
		}
	};
	template <class T, bool A>
	const T QuatT<T,A>::ZeroLen_Th;
	template <class T, bool A>
	const T QuatT<T,A>::Theta_Th;

	template <class T, bool A>
	inline std::ostream& operator << (std::ostream& os, const QuatT<T,A>& q) {
		os << "Quat: ";
		return q.asVec4().print(os);
	}
}
namespace std {
	template <class T, bool A>
	struct hash<frea::QuatT<T,A>> {
		using quat_t = frea::QuatT<T,A>;
		std::size_t operator()(const quat_t& q) const noexcept {
			using base_t = typename quat_t::base_t;
			return hash<base_t>()(static_cast<const base_t&>(q));
		}
	};
}
