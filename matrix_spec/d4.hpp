namespace frea {
	template <class V>
	struct MatT_spec<V,4,4> : MatT_dspec<V,4,4, MatT_spec<V,4,4>> {
		using this_t = MatT_spec;
		using base_t = MatT_dspec<V,4,4, MatT_spec<V,4,4>>;
		using base_t::base_t;
		using value_t = typename base_t::value_t;

		//! 透視変換行列
		static this_t PerspectiveFovLH(const Radian<value_t>& fov, const value_t& aspect, const value_t& nz, const value_t& fz);
		static this_t PerspectiveFovRH(const Radian<value_t>& fov, const value_t& aspect, const value_t& nz, const value_t& fz);
		static this_t _PerspectiveFov(const Radian<value_t>& fov, const value_t& aspect, const value_t& nz, const value_t& fz, const value_t& coeff);
	};
}
