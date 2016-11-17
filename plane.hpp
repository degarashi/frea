#pragma once
#include "matrix.hpp"

namespace frea {
	struct CantReachToPlane : std::runtime_error {
		CantReachToPlane(): std::runtime_error("CantReachToPlane") {}
	};

	template <class T, bool A>
	struct PlaneT :
		Data<T,4,A>,
		lubee::op::Mul<PlaneT<T,A>>,
		lubee::op::Ne<PlaneT<T,A>>
	{
		using op_m = lubee::op::Mul<PlaneT<T,A>>;
		using base_t = Data<T,4,A>;
		using base_t::base_t;
		using value_t = typename base_t::value_t;
		using vec_t = Vec_t<T,3,A>;
		using vec4_t = typename vec_t::template type_cn<4>;
		using mat3_t = Mat_t<T,3,3,A>;
		using mat4_t = typename mat3_t::template type_cn<4,4>;

		PlaneT() = default;
		constexpr PlaneT(const base_t& b) noexcept: base_t(b) {}
		// 違う内部形式からの変換(float->double など)
		template <class T2, bool A2>
		constexpr PlaneT(const PlaneT<T2,A2>& p) noexcept: base_t(p.asVec4()) {}
		constexpr PlaneT(const vec_t& orig, const value_t& dist) noexcept:
			base_t(orig.x, orig.y, orig.z, dist) {}

		//! 座標1つと法線から平面を作成
		constexpr static PlaneT FromPtDir(const vec_t& pos, const vec_t& dir) noexcept {
			return {dir, -dir.dot(pos)};
		}
		//! 座標3つから平面を作成
		constexpr static PlaneT FromPts(const vec_t& p0, const vec_t& p1, const vec_t& p2) noexcept {
			vec_t nml = (p1 - p0).cross(p2 - p0);
			nml.normalize();
			return {nml, -p0.dot(nml)};
		}
		//! 3つの平面が交差する地点を算出
		static auto ChokePoint(const PlaneT& p0, const PlaneT& p1, const PlaneT& p2) noexcept {
			struct {
				vec_t	pt;
				bool	cross;
			} ret;
			ret.cross = false;

			mat3_t m, mInv;
			m.template getRow<0>() = p0.getNormal();
			m.template getRow<1>() = p1.getNormal();
			m.template getRow<2>() = p2.getNormal();
			m.transpose();

			try {
				mInv = m.inversion();
			} catch(const NoInverseMatrix&) {
				return ret;
			}

			vec_t v(-p0.w, -p1.w, -p2.w);
			v *= mInv;
			ret.cross = true;
			ret.pt = v;
			return ret;
		}
		//! 2つの平面が交差する直線を算出
		static auto CrossLine(const PlaneT& p0, const PlaneT& p1) noexcept {
			struct {
				vec_t	pt, dir;
				bool	cross;
			} ret;
			ret.cross = false;

			const auto &nml0 = p0.getNormal(),
						&nml1 = p1.getNormal();
			vec_t nml = nml0 % nml1;
			if(std::abs(nml.len_sq()) < 1e-5)
				return ret;
			nml.normalize();

			mat3_t m, mInv;
			m.template getRow<0>() = nml0;
			m.template getRow<1>() = nml1;
			m.template getRow<2>() = nml;
			m.transpose();

			try {
				mInv = m.inversion();
			} catch(const NoInverseMatrix&) {
				return ret;
			}

			ret.cross = true;
			vec_t v(-p0.w, -p1.w, 0);
			v *= mInv;
			ret.pt = v;
			ret.dir = nml;
			return ret;
		}

		value_t dot(const vec_t& p) const noexcept {
			return asVec4().dot(p.template convertI<4,3>(1));
		}
		//! 平面を法線の方向へ指定距離、移動
		void move(const value_t& d) noexcept {
			this->w += d;
		}
		const vec_t& getNormal() const noexcept {
			return reinterpret_cast<const vec_t&>(*this);
		}
		template <class M>
		PlaneT operator * (const M& m) const {
			const auto& nml = getNormal();
			const vec_t tmp(nml * -this->w);
			return FromPtDir(
				(tmp.template convertI<4,3>(1)*m).template convert<3>(),
				((nml.template convert<4>()*m).normalization()).template convert<3>()
			);
		}
		using op_m::operator *;
		//! 引数の座標を垂直に平面へ下ろした地点の座標
		vec_t placeOnPlane(const vec_t& src) const noexcept {
			const value_t d = dot(src);
			return src - getNormal() * d;
		}
		//! 引数の座標をdir方向に向かって移動させ、平面と交差する座標
		value_t placeOnPlaneDirDist(const vec_t& dir, const vec_t& src) const {
			const value_t r = dir.dot(getNormal());
			if(std::abs(r) < 1e-6) {
				// dirと法線がほぼ平行なので幾ら頂点を移動させても平面へ辿りつけない
				throw CantReachToPlane();
			}
			return -dot(src) / r;
		}
		//! 平面上の一点
		vec_t getOrigin() const noexcept {
			return getNormal() * this->w;
		}
		//! Vector4に読み替え
		const vec4_t& asVec4() const noexcept {
			return reinterpret_cast<const vec4_t&>(*this);
		}
		bool operator == (const PlaneT& p) const noexcept {
			return asVec4() == p.asVec4();
		}
		//! 平面にてベクトルを反転
		vec_t flip(const vec_t& v) const noexcept {
			const auto& nml = getNormal();
			const value_t d = dot(v);
			return {v + nml*-d*2};
		}
		//! 平面との交差点を算出
		auto crosspoint(const vec_t& v0, const vec_t& v1) const noexcept {
			struct {
				vec_t point;
				bool cross;
			} res;
			// 線分が平面をまたぐか
			const value_t distf = dot(v0);
			const value_t distb = dot(v1);
			if(distf * distb >= 0) {
				res.cross = false;
				return res;
			}
			res.cross = true;
			const value_t div = std::abs(distf) + std::abs(distb);
			if(div < 1e-6) {
				// どちらの頂点ともほぼ面上にあるので交点はその中間ということにする
				res.point = (v0 + v1) /2;
			} else {
				const value_t ratio = std::abs(distf) / div;
				// 平面と線分の交点 -> tv
				res.point = (v1 - v0) * ratio + v0;
			}
			return res;
		}
	};

	using Plane = PlaneT<float, false>;
	using APlane = PlaneT<float, true>;
	using DPlane = PlaneT<double, false>;
	using ADPlane = PlaneT<double, true>;
}
namespace std {
	template <class T, bool A>
	struct hash<frea::PlaneT<T,A>> {
		using plane_t = frea::PlaneT<T,A>;
		std::size_t operator()(const plane_t& p) const noexcept {
			using base_t = typename plane_t::base_t;
			return hash<base_t>()(static_cast<const base_t&>(p));
		}
	};
}
