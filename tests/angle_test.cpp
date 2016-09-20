#include "test.hpp"
#include "../angle_func.hpp"

namespace frea {
	namespace test {
		template <class T>
		class Angle : public Random {
			protected:
				using value_t = T;
				using vec_t = Vec_t<value_t, 2, false>;
				using mat_t = Mat_t<value_t, 2, 2, false>;
				using deg_t = ::frea::Degree<T>;
				using rad_t = ::frea::Radian<T>;
				using Rd = decltype(std::declval<lubee::RandomMT>().template getUniformF<value_t>());
				Rd	_rd;

			public:
				Angle():
					_rd(mt().template getUniformF<value_t>({-1e3, 1e3}))
				{}
				T makeRF() {
					return _rd();
				}
				lubee::Range<value_t> makeRange() {
					return lubee::random::GenRange<value_t>(_rd);
				}
				rad_t makeRadian() {
					return rad_t(makeRF());
				}
				deg_t makeDegree() {
					return deg_t(makeRF());
				}
		};
		using Types = ToTestTypes_t<types::Float_t>;
		TYPED_TEST_CASE(Angle, Types);

		// ループ毎にsingleした値と独自にwhileで求めた値を比べる
		TYPED_TEST(Angle, Single) {
			constexpr auto onerotation = AngleInfo<Degree_t>::one_rotation<TypeParam>;
			auto angle = this->makeRF();
			typename TestFixture::deg_t degf(angle);
			constexpr int N_Iterations = 100;
			for(int i=0 ; i<N_Iterations ; i++) {
				auto tmp_degf = degf;
				{
					tmp_degf.single();

					auto val = angle;
					while(val >= onerotation)
						val -= onerotation;
					while(val < 0)
						val += onerotation;
					EXPECT_NEAR(tmp_degf.get(), val, lubee::ThresholdF<TypeParam>(0.5));
				}

				auto a = this->makeRF();
				angle += a;
				degf.set(degf.get() + a);
			}
		}
		TYPED_TEST(Angle, Single_Semicircle) {
			auto a0 = this->makeDegree(),
				 a1 = a0;
			a0.single();
			a0.semicircle();
			a1.semicircle();
			a1.single();
			ASSERT_LT(AngleDiff(a0, a1).get(), lubee::ThresholdF<TypeParam>(0.6));
		}
		TYPED_TEST(Angle, Range) {
			USING(deg_t);
			deg_t degf(this->makeRF());
			auto r = this->makeRange();
			auto val = Saturate(degf.get(), r.from, r.to);
			degf.range({Degree<TypeParam>(r.from), Degree<TypeParam>(r.to)});
			EXPECT_EQ(degf.get(), val);
		
			// rangeValueテスト
			degf.set(this->makeRF());
			r = this->makeRange();
			val = Saturate(degf.get(), r.from, r.to);
			degf.rangeValue(r);
			EXPECT_EQ(degf.get(), val);
		}
		TYPED_TEST(Angle, Degree_Radian_Degree) {
			// Degree -> Radian -> Degree で値の比較
			const Degree<TypeParam> degf(this->makeRF());
			const auto radf = degf.template convert<Radian_t>();
			const auto degf2 = radf.template convert<Degree_t>();
			EXPECT_NEAR(degf.get(), degf2.get(), lubee::ThresholdF<TypeParam>(0.6));
		}
		TYPED_TEST(Angle, Degree_Radian_Convert) {
			const auto deg = this->makeRF();
			const Degree<TypeParam> degf(deg);
			const Radian<TypeParam> radf(degf);
			const auto r0 = degf.get() / AngleInfo<Degree_t>::one_rotation<TypeParam>;
			const auto r1 = radf.get() / AngleInfo<Radian_t>::one_rotation<TypeParam>;
			EXPECT_NEAR(r0, r1, lubee::ThresholdF<TypeParam>(0.6));
		}
		TYPED_TEST(Angle, Arithmetic) {
			USING(deg_t);
			deg_t	degf(this->makeRF()),
					degf2(this->makeRF());
			EXPECT_EQ((degf + degf2).get(), degf.get() + degf2.get());
			EXPECT_EQ((degf - degf2).get(), degf.get() - degf2.get());
			EXPECT_EQ((degf * 2).get(), degf.get() * 2);
			EXPECT_EQ((degf / 2).get(), degf.get() / 2);
		
			auto val = static_cast<TypeParam>(degf) + static_cast<TypeParam>(degf2);
			EXPECT_EQ((degf += degf2).get(), val);
			val = static_cast<TypeParam>(degf) - static_cast<TypeParam>(degf2);
			EXPECT_EQ((degf -= degf2).get(), val);
			val = static_cast<TypeParam>(degf) * 2;
			EXPECT_EQ((degf *= 2).get(), val);
			val = static_cast<TypeParam>(degf) / 2;
			EXPECT_EQ((degf /= 2).get(), val);
		}
		TYPED_TEST(Angle, Lerp) {
			USING(deg_t);
			deg_t	 ang0(this->makeRF()),
					ang1(this->makeRF()),
					tmp, tmp2;
			ang0.single();
			ang1.single();
			// Lerp係数が0ならang0と等しい
			tmp = AngleLerp(ang0, ang1, 0.0);
			tmp.single();
			constexpr auto Th = lubee::ThresholdF<TypeParam>(0.5);
			EXPECT_NEAR(ang0.get(), tmp.get(), Th);
			// Lerp係数が1ならang1と等しい
			tmp = AngleLerp(ang0, ang1, 1.0);
			tmp.single();
			EXPECT_NEAR(ang1.get(), tmp.get(), Th);
			// Lerp係数を0.5で2回かければ0.75でやったのと(ほぼ)同じになる
			tmp = AngleLerp(ang0, ang1, .5);
			tmp = AngleLerp(tmp, ang1, .5);
			tmp2 = AngleLerp(ang0, ang1, .75);
			tmp.single(); tmp2.single();
			EXPECT_NEAR(tmp.get(), tmp2.get(), Th);
		}
		TYPED_TEST(Angle, Move) {
			USING(value_t);
			USING(deg_t);
			deg_t	ang0(this->makeRF()),
					ang1;
			constexpr auto oa = deg_t::OneRotationAng;
			ang0.single();
			value_t diff;
			do {
				ang1.set(this->makeRF());
				ang1.single();
				diff = AngleLerpValueDiff(ang0.get(), ang1.get(), oa);
			} while(diff < 1e-2);
			deg_t tmp = ang0;
			// 角度の差分を3で割ったら3回処理した時点でang1とほぼ等しくなる
			const auto diff3 = deg_t(std::abs(diff / 3));
			auto fnDiff = [&tmp, &ang1, oa](){
				return std::abs(AngleLerpValueDiff(tmp.get(), ang1.get(), oa));
			};
			tmp = AngleMove(tmp, ang1, diff3);
			constexpr auto Th = lubee::ThresholdF<TypeParam>(0.6);
			EXPECT_FALSE(fnDiff() < Th);
			tmp = AngleMove(tmp, ang1, diff3);
			EXPECT_FALSE(fnDiff() < Th);
			tmp = AngleMove(tmp, ang1, diff3);
			EXPECT_TRUE(fnDiff() < Th);
		}
		TYPED_TEST(Angle, AngleValue) {
			USING(deg_t);
			USING(vec_t);
			USING(mat_t);
			deg_t ang(this->makeRF());
			// VectorFromAngleで計算したものと行列で計算したものとはほぼ一致する
			const vec_t v0 = VectorFromAngle(ang),
						v1 = vec_t(0,1) * mat_t::Rotation(ang);
			constexpr auto Th = lubee::ThresholdF<TypeParam>(0.6);
			EXPECT_NEAR(v0.x, v1.x, Th);
			EXPECT_NEAR(v0.y, v1.y, Th);

			// AngleValueにかければ元の角度が出る
			deg_t ang1 = AngleValue(v0);
			ang.single();
			ang1.single();
			EXPECT_LT(AngleDiff(ang, ang1).get(), lubee::ThresholdF<TypeParam>(0.9));
		}
		TYPED_TEST(Angle, Serialization) {
			USING(deg_t);
			USING(rad_t);
			lubee::CheckSerialization(deg_t(this->makeRF()));
			lubee::CheckSerialization(rad_t(this->makeRF()));
		}
	}
}
