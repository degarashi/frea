namespace frea {
	template <class V>
	struct MatT_spec<V,3,3> : MatT_dspec<V,3,3, MatT_spec<V,3,3>> {
		using this_t = MatT_spec;
		using base_t = MatT_dspec<V,3,3, MatT_spec<V,3,3>>;
		using base_t::base_t;
		using value_t = typename base_t::value_t;
		using vec_t = typename base_t::vec_t;

		//! X軸周りの回転
		static this_t RotationX(RadF ang);
		//! Y軸周りの回転
		static this_t RotationY(RadF ang);
		//! Z軸周りの回転
		static this_t RotationZ(RadF ang);
		//! 任意軸周りの回転
		static this_t RotationAxis(const vec_t& axis, const Radian<value_t>& ang);
		static this_t LookAtLH(const vec_t& pos, const vec_t& at, const vec_t& up);
		static this_t LookDirLH(const vec_t& pos, const vec_t& at, const vec_t& up);
		static this_t LookAtRH(const vec_t& pos, const vec_t& at, const vec_t& up);
		static this_t LookDirRH(const vec_t& pos, const vec_t& at, const vec_t& up);
		this_t transposition() const { return this->asInternal().transposition(); }
		this_t transpose() { return this->asInternal().transpose(); }
	};
}
