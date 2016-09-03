#include "test.hpp"

namespace frea {
	namespace test {
		template <class T>
		using MatrixD_4 = RMatrix<T>;
		using TypesD_4 = ToTestTypes_t<types::SMatrixRange_t<types::FReg_t, 4,5>>;
		TYPED_TEST_CASE(MatrixD_4, TypesD_4);

		TYPED_TEST(MatrixD_4, LookAt) {
			USING(mat_t);
			USING(value_t);
			using vec3_t = typename mat_t::vec3_t;

			constexpr auto range = Range<value_t>{-1e2, 1e2};
			const auto mtf = this->mt().template getUniformF<value_t>(range);
			// {pos, at}
			const auto pos = MakePos<2, vec3_t>(mtf, 1e-2);
			const vec3_t dir = (pos[1]-pos[0]).normalization();
			vec3_t up;
			do {
				up = this->makeDir();
			} while(std::abs(up.dot(dir)) < 0.9);
			up.normalize();
			const auto m0 = mat_t::LookAt(pos[0], pos[1], up),
						m1 = mat_t::LookDir(pos[0], dir, up);
			constexpr auto Th = ThresholdF<value_t>(0.7);
			// ビュー行列で変換した座標点(at)はZ軸上にある(X==0, Y==0, Z=distance(v - pos))
			const vec3_t v_at[2] = {
				(pos[1].template convertI<4,3>(1)*m0).template convert<3>(),
				(pos[1].template convertI<4,3>(1)*m1).template convert<3>()
			};
			for(auto& vp : v_at) {
				EXPECT_NEAR(vp.x, 0, Th);
				EXPECT_NEAR(vp.y, 0, Th);
				EXPECT_NEAR(vp.z, pos[1].distance(pos[0]), Th);
			}
			// 任意のビュー行列を2通りの方法で作っても結果は同じ
			EXPECT_LE(AbsMax(mat_t(m1-m0)), ThresholdF<value_t>(0.3));
			// AtとPosが同一座標の場合はNoValidAxis例外を送出
			EXPECT_THROW(mat_t::LookAt(pos[0], pos[0], up), NoValidAxis);
			// upとdirの内積が1または-1の場合、NoValidAxis例外を送出
			EXPECT_THROW(mat_t::LookDir(pos[0], up, up), NoValidAxis);
			EXPECT_THROW(mat_t::LookDir(pos[0], -up, up), NoValidAxis);
		}
		TYPED_TEST(MatrixD_4, PerspectiveFov) {
			USING(value_t);
			USING(rad_t);
			USING(mat_t);
			USING(vec_t);
			using vec3_t = typename mat_t::vec3_t;

			const auto mtf = this->mt().template getUniformF<value_t>();
			const rad_t fov = Degree<value_t>(mtf({10,150}));
			const value_t nz = mtf({1e-3, 1e3}),
						fz = nz+mtf({1e-3, 1e3}),
						aspect = mtf({1e-1, 1e1});
			const auto m = mat_t::PerspectiveFov(fov, aspect, nz, fz);
			const value_t xr = mtf({-1,1}),
						yr = mtf({-1,1}),
						zr = mtf({0,1});
			const auto ft2 = std::tan(fov.get()/2);
			const vec_t v0(ft2*aspect * xr,
							ft2 * yr,
							(fz-nz)*zr+nz,
							1),
						v1 = v0 * m;

			// (X,Y,Z,1) -> (dir.x*aspect*tan(fov/2), dir.y*tan(fov/2), [0->W], W)
			const vec3_t v2 = v1.asVec3Coord();
			EXPECT_TRUE(IsInRange<value_t>(v2.x, -1, 1));
			EXPECT_TRUE(IsInRange<value_t>(v2.y, -1, 1));
			EXPECT_TRUE(IsInRange<value_t>(v2.z, 0, 1));
		}
	}
}
