#pragma once
#include "matrix.hpp"

namespace frea {
	template <class T, bool A>
	struct PlaneT : Data<T,4,A> {
		using base_t = Data<T,4,A>;
		using base_t::base_t;
		using value_t = typename base_t::value_t;
		using vec_t = Vec_t<T,3,A>;
		using mat3_t = MatT<vec_t, 3>;
		using mat4_t = MatT<typename vec_t::template type_cn<4>, 4>;

		PlaneT() = default;
		constexpr PlaneT(const base_t& b):
			base_t(b) {}
		// 違う要素Quatからの変換
		template <class T2, bool A2>
		constexpr PlaneT(const PlaneT<T2,A2>& p):
			base_t(p.asVec4()) {}
		constexpr PlaneT(const vec_t& orig, const value_t& dist):
			base_t(orig.x, orig.y, orig.z, dist) {}

		constexpr static PlaneT FromPtDir(const vec_t& pos, const vec_t& dir) {
			return {dir, -dir.dot(pos)};
		}
		constexpr static PlaneT FromPts(const vec_t& p0, const vec_t& p1, const vec_t& p2) {
			vec_t nml = (p1 - p0).cross(p2 - p0);
			nml.normalize();
			return {nml, -p0.dot(nml)};
		}
		constexpr static std::tuple<vec_t,bool> ChokePoint(const PlaneT& p0, const PlaneT& p1, const PlaneT& p2) {
			mat3_t m, mInv;
			m.getRow(0) = p0.getNormal();
			m.getRow(1) = p1.getNormal();
			m.getRow(2) = p2.getNormal();
			m.transpose();

			if(!m.inversion(mInv))
				return std::make_tuple(vec_t(), false);

			vec_t v(-p0.d, -p1.d, -p2.d);
			v *= mInv;
			return std::make_tuple(v, true);
		}
		constexpr static std::tuple<vec_t,vec_t,bool> CrossLine(const PlaneT& p0, const PlaneT& p1) {
			const auto &nml0 = p0.getNormal(),
						&nml1 = p1.getNormal();
			vec_t nml = nml0 % nml1;
			if(std::fabs(nml.len_sq()) < 1e-5)
				return std::make_tuple(vec_t(),vec_t(), false);
			nml.normalize();

			mat3_t m, mInv;
			m.getRow(0) = nml0;
			m.getRow(1) = nml1;
			m.getRow(2) = nml;
			m.transpose();

			if(!m.inversion(mInv))
				return std::make_tuple(vec_t(),vec_t(),false);

			vec_t v(-p0.d, -p1.d, 0);
			v *= mInv;
			return std::make_tuple(v, nml, true);
		}

		value_t dot(const vec_t& p) const {
			const auto v = asVec4() * p.template convert<4>();
			return v.sumUp() + this->w;
		}
		void move(const value_t& d) {
			this->w += d;
		}
		const vec_t& getNormal() const {
			return static_cast<const vec_t&>(*this);
		}
		PlaneT operator * (const mat4_t& m) const {
			const auto& nml = getNormal();
			const vec_t tmp(nml * -this->w);
			return FromPtDir(tmp.convertI<4>(1)*m, (nml.template convert<4>()*m).normalization());
		}
		PlaneT& operator *= (const mat4_t& m) {
			return *this = *this * m;
		}
		vec_t placeOnPlane(const vec_t& src, const value_t& offset) const {
			const value_t d = dot(src);
			return src + getNormal() * (-this->w+offset);
		}
		value_t placeOnPlaneDirDist(const vec_t& dir, const vec_t& src) const {
			const value_t r = dir.dot(getNormal());
			if(std::fabs(r) < 1e-6)
				throw 0;
			const auto v = src - getOrigin();
			const value_t act = dot(v);
			return act / r;
		}
		vec_t getOrigin() const {
			return getNormal() * -this->w;
		}
		const vec_t& asVec4() const noexcept {
			return static_cast<const vec_t&>(*this);
		}
		bool operator == (const PlaneT& p) const {
			return static_cast<const base_t&>(*this)
					== static_cast<const base_t&>(p);
		}
		bool operator != (const PlaneT& p) const {
			return !(*this == p);
		}
	};
}
