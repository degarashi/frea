#include "test.hpp"
#include "../plane.hpp"
#include "../random/plane.hpp"

namespace frea {
	namespace test {
		template <class T>
		class Plane : public Random {
			public:
				using value_t = typename std::tuple_element_t<0,T>;
				constexpr static bool align = std::tuple_element_t<1,T>::value;
				using plane_t = PlaneT<value_t, align>;
				using vec_t = typename plane_t::vec_t;
			private:
				constexpr static Range<value_t> DefaultRange{-1e2, 1e2};
				using RD = decltype(std::declval<Random>().mt().template getUniformF<value_t>(DefaultRange));
				RD	_rd;
			public:
				Plane():
					_rd(mt().template getUniformF<value_t>(DefaultRange))
				{}
				const RD& getMtf() const {
					return _rd;
				}
				auto makePlane() const {
					return random::GenPlane<plane_t>(_rd);
				}
				auto makeVec3() const {
					return random::GenVec<vec_t>(_rd);
				}
				auto makeDir() const {
					return random::GenVecUnit<vec_t>(_rd);
				}
				auto makeRF() const {
					return _rd();
				}
				auto makeAxis() const {
					struct {
						vec_t x,y,z;
					} ret = {{1,0,0}, {0,1,0}, {0,0,1}};
					const auto q = random::GenQuat<QuatT<value_t,false>>(_rd);
					ret.x *= q;
					ret.y *= q;
					ret.z *= q;
					return ret;
				}
				auto makeOrthogonalMat() const {
					const auto axis = makeAxis();
					Mat_t<value_t, 4, 4, false> m;
					m.template setColumn<0>(axis.x.template convertI<4,3>(makeRF()));
					m.template setColumn<1>(axis.y.template convertI<4,3>(makeRF()));
					m.template setColumn<2>(axis.z.template convertI<4,3>(makeRF()));
					m.template setColumn<3>({0,0,0,1});
					return m;
				}
				// 位置が重複しないランダム頂点配列
				auto makeVec3s(const int n) const {
					std::vector<vec_t> ret(n);
					ret[0] = makeVec3();
					for(int i=1 ; i<n ; i++) {
						for(;;) {
							ret[i] = makeVec3();
							bool b = true;
							for(int j=0 ; j<i ; j++) {
								if(ret[j].dist_sq(ret[i]) < 1e-6) {
									b = false;
									break;
								}
							}
							if(b)
								break;
						}
					}
					return ret;
				}
		};
		template <class T>
		constexpr Range<typename Plane<T>::value_t> Plane<T>::DefaultRange;

		TYPED_TEST_CASE(Plane, types::PTypes);

		TYPED_TEST(Plane, Compare) {
			const auto p0 = this->makePlane();
			auto p1 = p0;

			// 自身との==はTrue
			EXPECT_EQ(p0, p0);
			// コピーした変数との==もTrue
			EXPECT_EQ(p0, p1);
			// !=は==の反対
			EXPECT_FALSE(p0 != p1);

			// 値を1ついじった場合は==がFalse
			auto& mt = this->mt();
			const int idx = mt.template getUniform<int>({0,3});
			p1[idx] += 1;
			EXPECT_NE(p0, p1);
		}
		TYPED_TEST(Plane, Generate) {
			USING(plane_t);
			USING(value_t);

			constexpr auto Th = ThresholdF<value_t>(0.6);
			const auto vtx = this->makeVec3s(3);
			{
				auto p = plane_t::FromPts(vtx[0], vtx[1], vtx[2]);
				// 3つの頂点は平面上にある
				for(int i=0 ; i<3 ; i++)
					ASSERT_LE(std::abs(p.dot(vtx[i])), Th);

				// 平面を少しずらせば頂点群はその距離だけ平面から離れる
				const auto dist = this->makeRF();
				p.move(dist);
				for(int i=0 ; i<3 ; i++)
					ASSERT_LE(p.dot(vtx[i])-dist, Th);
			}
			{
				const auto dir = this->makeDir();
				const auto p = plane_t::FromPtDir(vtx[0], dir);
				// vtx[0]は平面上
				ASSERT_LE(std::abs(p.dot(vtx[0])), Th);
			}
		}
		TYPED_TEST(Plane, Transform) {
			USING(value_t);
			// 頂点と平面を行列(直交行列)で変換しても相対的な位置関係は変わらない
			auto p = this->makePlane();
			auto v = this->makeVec3();
			const auto d0 = p.dot(v);
			const auto m = this->makeOrthogonalMat();

			p *= m;
			v = (v.template convertI<4,3>(1) * m).template convert<3>();
			const auto d1 = p.dot(v);

			constexpr auto Th = ThresholdF<value_t>(0.6);
			ASSERT_LE(std::abs(d1-d0), Th);
		}
		TYPED_TEST(Plane, Place) {
			USING(vec_t);
			USING(value_t);
			auto p = this->makePlane();
			const auto v0 = this->makeVec3();
			{
				const auto v1 = p.placeOnPlane(v0);
				// 平面との距離はゼロである
				ASSERT_LE(std::abs(p.dot(v1)), ThresholdF<value_t>(0.5));
			}
			{
				p.x = 1;
				p.y = 0;
				p.z = 0;
				p.w = 1;
				vec_t dir;
				do {
					// 平面に平行でない方向ベクトルが出るまでランダム生成
					dir = this->makeDir();
					dir = vec_t(-1,0,0).normalization();
				} while(std::abs(dir.dot(p.getNormal())) < 0.25);
				const auto dist = p.placeOnPlaneDirDist(dir, v0);
				const auto v1 = v0 + dir*dist;
				ASSERT_LE(std::abs(p.dot(v1)), ThresholdF<value_t>(0.5));
			}
		}
		TYPED_TEST(Plane, CrossLine) {
			USING(plane_t);
			USING(value_t);
			// 互いに平行でない平面を2つ生成
			const plane_t p0 = this->makePlane();
			plane_t p1;
			do {
				p1 = this->makePlane();
			} while(std::abs(p0.getNormal().dot(p1.getNormal())) < 0.25);

			// 交差する部分は直線になるので、その上で点を移動させつつ双方の面上にあるか確認
			const auto res = plane_t::CrossLine(p0, p1);
			ASSERT_TRUE(res.cross);

			// チェックする終端を適当に決める
			constexpr int NDiv = 8;
			for(int i=0 ; i<NDiv ; i++) {
				const auto pos = res.pt + res.dir * value_t(i)/NDiv;
				constexpr auto Th = ThresholdF<value_t>(0.8);
				ASSERT_LE(std::abs(p0.dot(pos)), Th);
				ASSERT_LE(std::abs(p1.dot(pos)), Th);
			}
		}
		TYPED_TEST(Plane, ChokePoint) {
			USING(plane_t);
			USING(value_t);
			USING(vec_t);
			// 互いに平行でない平面を3つ生成
			plane_t p[3];
			p[0] = this->makePlane();
			for(int i=1 ; i<3 ; i++) {
				for(;;) {
					p[i] = this->makePlane();
					bool b = true;
					for(int j=0 ; j<i ; j++) {
						if(std::abs(p[i].getNormal().dot(p[j].getNormal())) < 0.5) {
							b = false;
							break;
						}
					}
					if(b)
						break;
				}
			}

			const auto res = plane_t::ChokePoint(p[0], p[1], p[2]);
			if(res.cross) {
				// 交差する部分は点になるので、その点が全ての面上にあるか確認
				const auto Th = res.pt.length() / 4096;
				for(auto& pt : p)
					ASSERT_LE(std::abs(pt.dot(res.pt)), Th);
			} else {
				// 3つの平面が一箇所で交わる点がない -> それぞれ2つを選んだ時の交点（直線）が平行になる
				vec_t dir[3];
				int ndir = 0;
				for(int i=0 ; i<3 ; i++) {
					const auto res = plane_t::CrossLine(p[i], p[(i+1)%3]);
					if(res.cross)
						dir[ndir++] = res.dir;
				}
				constexpr auto Th = ThresholdF<value_t>(0.6);
				// 交わる平面の組み合わせが少なくとも2つはある
				ASSERT_GE(ndir, 2);
				for(int i=0 ; i<ndir-1 ; i++) {
					ASSERT_LE(std::abs(dir[i].dot(dir[i+1]))-1, Th);
				}
			}
		}
		TYPED_TEST(Plane, Serialization) {
			CheckSerialization(this->makePlane());
		}
	}
}
