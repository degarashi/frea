namespace frea {
	template <class V>
	struct MatT_spec<V,3,3> : MatT_dspec<V,3,3, MatT_spec<V,3,3>> {
		using this_t = MatT_spec;
		using base_t = MatT_dspec<V,3,3, MatT_spec<V,3,3>>;
		using base_t::base_t;
		using value_t = typename base_t::value_t;
		using vec_t = typename base_t::vec_t;
		using rad_t = Radian<value_t>;

		//! X軸周りの回転
		static this_t RotationX(const rad_t& ang) noexcept {
			const value_t	C = std::cos(ang.get()),
							S = std::sin(ang.get());
			return {
				1,	0,	0,
				0,	C,	-S,
				0,	S,	C,
			};
		}
		//! Y軸周りの回転
		static this_t RotationY(const rad_t& ang) noexcept {
			const value_t	C = std::cos(ang.get()),
							S = std::sin(ang.get());
			return {
				C,	0,	S,
				0,	1,	0,
				-S,	0,	C
			};
		}
		//! Z軸周りの回転
		static this_t RotationZ(const rad_t& ang) noexcept {
			const value_t	C = std::cos(ang.get()),
							S = std::sin(ang.get());
			return {
				C,	-S,	0,
				S,	C,	0,
				0,	0,	1
			};
		}
		//! 任意軸周りの回転
		static this_t RotationAxis(const vec_t& axis, const rad_t& ang) noexcept {
			const value_t	C = std::cos(ang.get()),
							S = std::sin(ang.get()),
							RC = 1-C;
			const value_t M00 = C + Square(axis.x) * RC,
						M01 = axis.x * axis.y * RC + axis.z*S,
						M02 = axis.x * axis.z * RC - axis.y*S,
						M10 = axis.x * axis.y * RC - axis.z*S,
						M11 = C + Square(axis.y) * RC,
						M12 = axis.y * axis.z * RC + axis.x*S,
						M20 = axis.x * axis.z * RC + axis.y*S,
						M21 = axis.y * axis.z * RC - axis.x*S,
						M22 = C + Square(axis.z) * RC;
			return {
				M00, M10, M20,
				M01, M11, M21,
				M02, M12, M22
			};
		}
		static this_t LookAtLH(const vec_t& pos, const vec_t& at, const vec_t& up);
		static this_t LookDirLH(const vec_t& pos, const vec_t& at, const vec_t& up);
		static this_t LookAtRH(const vec_t& pos, const vec_t& at, const vec_t& up);
		static this_t LookDirRH(const vec_t& pos, const vec_t& at, const vec_t& up);
		this_t transposition() const { return this->asInternal().transposition(); }
		void transpose() { *this = transposition(); }
	};
}