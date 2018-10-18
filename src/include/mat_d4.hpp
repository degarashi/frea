namespace frea {
	template <class V>
	struct MatT_spec<V,4,4> : MatT_dspec<V,4,4, MatT_spec<V,4,4>> {
		using this_t = MatT_spec;
		using base_t = MatT_dspec<V,4,4, MatT_spec<V,4,4>>;
		using base_t::base_t;
		using value_t = typename base_t::value_t;
		using vec_t = typename base_t::vec_t;
		using vec3_t = typename vec_t::template type_cn<3>;

		//! 位置と注視点からビュー行列生成
		static this_t LookAt(const vec3_t& pos, const vec3_t& at, const vec3_t& up) {
			return LookDir(pos, (at-pos).normalization(), up);
		}
		//! 位置と注視方向からビュー行列生成
		static this_t LookDir(const vec3_t& pos, const vec3_t& dir, const vec3_t& up) {
			vec3_t xA(up % dir),
				   yA;
			if(xA.normalize() < lubee::Threshold<typename vec3_t::value_t>(0.3, 1) ||
				xA.isNaN())
				throw NoValidAxis();
			yA = dir % xA;

			return this_t(
				xA.x, yA.x, dir.x, 0,
				xA.y, yA.y, dir.y, 0,
				xA.z, yA.z, dir.z, 0,
				-pos.dot(xA), -pos.dot(yA), -pos.dot(dir), 1
			);
		}
		using Range_t = lubee::Range<value_t>;
		using Radian_t = Radian<value_t>;
		//! 透視変換行列
		static this_t PerspectiveFov(const Radian_t& fov, const value_t& aspect, const Range_t& z) {
			return _PerspectiveFov(fov, aspect, z, 1);
		}
		static this_t _PerspectiveFov(const Radian_t& fov, const value_t& aspect, const Range_t& z, const value_t& coeff) {
			if(fov.get() <= 0 ||
				fov.get() >= Radian_t::HalfRotationRange.to)
				throw InvalidFov();
			value_t h = 1.0f / std::tan(fov.get()/2),
					w = h / aspect,
					f0 = z.to/(z.to-z.from),
					f1 = -z.from*z.to/(z.to-z.from);
			return this_t(
				w,	0,		0,			0,
				0,	h,		0,			0,
				0,	0,		f0,			coeff,
				0,	0,		f1*coeff,	0
			);
		}
	};
}
