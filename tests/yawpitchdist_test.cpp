#include "test.hpp"
#include "../yawpitchdist.hpp"
#include "../serialization/yawpitchdist.hpp"

namespace frea {
	namespace test {
		template <class T>
		struct YawPitchDist : Random {
			using value_t = T;
		};
		using Types = ToTestTypes_t<types::Float_t>;
		TYPED_TEST_CASE(YawPitchDist, Types);

		TYPED_TEST(YawPitchDist, FromPos_Random) {
			USING(value_t);
			const auto mtf = this->mt().template getUniformF<value_t>({-1e2, 1e2});
			// YawPitchDistで計算した値を再度合成して結果を比較
			using V3 = Vec_t<value_t, 3, true>;
			V3 nzv;
			do {
				nzv = random::GenVec<V3>(mtf);
			} while(nzv.length() < 1e-2);
			const auto ypd = ::frea::YawPitchDist<value_t>::FromPos(nzv);
			const auto res = ypd.toOffsetRot();
			const auto nzv2 = -res.rot.getDir() * ypd.distance;

			constexpr auto Th = lubee::ThresholdF<value_t>(0.8);
			EXPECT_LT(AbsMax(V3(res.pos - nzv)), Th);
			// rotation + distanceから復元した元座標ベクトル
			EXPECT_LT(AbsMax(V3(nzv2 - nzv)), Th);
		}

		// XZ -> YawAngle
		struct Yaw_Angle {
			Vec3	dir;
			DegF	angle;
		};
		const Yaw_Angle c_yawAngle[] = {
			{{0,0,1}, DegF(0)},
			{{1,0,1}, DegF(45)},
			{{1,0,0}, DegF(90)},
			{{1,0,-1}, DegF(135)},
			{{0,0,-1}, DegF(180)},
			{{-1,0,-1}, DegF(225)},
			{{-1,0,0}, DegF(270)},
			{{-1,0,1}, DegF(315)}
		};
		TYPED_TEST(YawPitchDist, FromPos_Fixed) {
			USING(value_t);
			using deg_t = Degree<value_t>;
			const auto mtf = this->mt().template getUniformF<value_t>({-1,1});
			// 予め答えが用意されたパターンと照らし合わせ
			for(auto& ya : c_yawAngle) {
				auto pos = ya.dir;
				pos.y = mtf();	// Pitchはランダムな値
				const auto ypd = ::frea::YawPitchDist<value_t>::FromPos(pos);
				EXPECT_NEAR(deg_t(ya.angle).get(), deg_t(ypd.yaw).get(), lubee::ThresholdF<value_t>(0.8));
			}
		}
		TYPED_TEST(YawPitchDist, Serialization) {
			USING(value_t);
			using rad_t = Radian<value_t>;
			const auto mtf = this->mt().template getUniformF<value_t>();
			::frea::YawPitchDist<value_t> ypd;
			ypd.yaw = random::GenAngle<rad_t>(mtf);
			ypd.pitch = random::GenHalfAngle<rad_t>(mtf);
			ypd.distance = mtf({1e-2, 1e2});
			lubee::CheckSerialization(ypd);
		}
	}
}
