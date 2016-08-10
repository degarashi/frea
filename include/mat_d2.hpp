namespace frea {
	template <class V>
	struct MatT_spec<V,2,2> : MatT_dspec<V,2,2, MatT_spec<V,2,2>> {
		using base_t = MatT_dspec<V,2,2, MatT_spec<V,2,2>>;
		using base_t::base_t;
		using this_t = MatT_spec;
		using value_t = typename V::value_t;
		using vec_t = typename base_t::vec_t;

		static this_t Rotation(const Radian<value_t>& ang) {
			const value_t S = std::sin(ang.get()),
							C = std::cos(ang.get());
			return {
				C, S,
				-S, C
			};
		}
	};
}
