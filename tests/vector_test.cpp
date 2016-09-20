#include "test.hpp"
#include "../compiler_macro.hpp"

namespace frea {
	namespace test {
		// Int, Float両方のテストケース
		template <class T>
		using Vector = RVector<T>;
		using Types = ToTestTypes_t<types::VectorRange_t<types::Value_t>>;
		TYPED_TEST_CASE(Vector, Types);

		namespace {
			template <class T>
			constexpr lubee::Range<T> DefaultRange{-1e3, 1e3};
		}
		TYPED_TEST(Vector, Wrap_Equality) {
			USING(value_t);
			constexpr auto range = DefaultRange<value_t>;
			const auto v0 = this->makeRVec(range),
						v1 = this->makeRVec(range);
			const auto w0 = v0.asInternal(),
						w1 = v1.asInternal();
			const auto mtf = this->mt().template getUniformF<value_t>(range);
			auto s0 = mtf(),
				 s1 = mtf();
			if(s0 > s1)
				std::swap(s0, s1);

			// (vec +-*/ vec) == (wrap +-*/ wrap)
			EXPECT_EQ(v0+v1, w0+w1);
			EXPECT_EQ(v0-v1, w0-w1);
			EXPECT_EQ(v0*v1, w0*w1);
			// ゼロ除算避け
			constexpr auto Th_z = lubee::Threshold<value_t>(0.4, 1);
			if(!HasZero(v1.m, Th_z)) {
				EXPECT_EQ(v0/v1, w0/w1);
			}
			// (vec +-*/ s) == (wrap +-*/ s)
			EXPECT_EQ(v0+s0, w0+s0);
			EXPECT_EQ(v0-s0, w0-s0);
			EXPECT_EQ(v0*s0, w0*s0);
			if(std::abs(s0) >= Th_z) {
				EXPECT_EQ(v0/s0, w0/s0);
			}
			// (vec == vec) == (wrap == wrap)
			EXPECT_EQ(v0==v0, w0==w0);
			EXPECT_EQ(v0!=v0, w0!=w0);
			// dot
			EXPECT_EQ(v0.dot(v1), w0.dot(w1));
			// average
			EXPECT_EQ(v0.average(), w0.average());
			// saturation
			EXPECT_EQ(v0.saturation(s0, s1), w0.saturation(s0, s1));
			// getMin
			EXPECT_EQ(v0.getMin(v1), w0.getMin(w1));
			{
				// selectMin
				auto tv0 = v0;
				auto tw0 = w0;
				tv0.selectMin(v1);
				tw0.selectMin(w1);
				EXPECT_EQ(tv0, tw0);
			}
			// getMax
			EXPECT_EQ(v0.getMax(v1), w0.getMax(w1));
			{
				// selectMax
				auto tv0 = v0;
				auto tw0 = w0;
				tv0.selectMax(v1);
				tw0.selectMax(w1);
				EXPECT_EQ(tv0, tw0);
			}
			// op -
			EXPECT_EQ(-v0, -w0);
			// getMinValue
			EXPECT_EQ(v0.getMinValue(), w0.getMinValue());
			// getMaxValue
			EXPECT_EQ(v0.getMaxValue(), w0.getMaxValue());
			// absolute
			EXPECT_EQ(v0.absolute(), w0.absolute());
			{
				// isZero
				const auto th = this->mt().template getUniform<value_t>(range*2);
				EXPECT_EQ(v0.isZero(th), w0.isZero(th));
			}
		}
		TYPED_TEST(Vector, Iterator) {
			USING(value_t);
			// イテレータを介して値を加算した場合とインデックスを介した場合で比べる　
			constexpr auto R = DefaultRange<value_t>;
			auto v0 = this->makeRVec(R);
			auto v1 = v0;
			const auto value = this->mt().template getUniform<value_t>(R);
			for(auto itr=v1.begin() ; itr!=v1.end() ; ++itr)
				*itr += value;
			constexpr int N = decltype(v0)::size;
			for(int i=0 ; i<N ; i++)
				v0[i] += value;

			for(int i=0 ; i<N ; i++)
				ASSERT_FLOAT_EQ(v0[i], v1[i]);

			auto itr0 = v0.cbegin();
			auto itr1 = v1.cbegin();
			for(;;) {
				bool b0 = itr0==v0.cend(),
					 b1 = itr1==v1.cend();
				ASSERT_EQ(b0, b1);
				if(b0)
					break;
				ASSERT_FLOAT_EQ(*itr0, *itr1);
				++itr0;
				++itr1;
			}
		}
		TYPED_TEST(Vector, Minus) {
			USING(value_t);
			USING(array_t);
			USING(vec_t);

			// operator -
			constexpr auto range = DefaultRange<value_t>;
			auto vec = this->makeRVec(range);
			// 地道に要素の符号反転した結果と比べる
			array_t ar(vec);
			for(auto& a : ar.m)
				a = -a;
			vec = -vec;
			constexpr auto Th = lubee::Threshold<value_t>(0.1, 0);
			ASSERT_LE(AbsMax(ar - vec), Th);			// Check: Vectorクラス
			ASSERT_LE(AbsMax(ar - vec_t(vec+0)), Th);	// Check: Wrapクラス
		}
		namespace {
			template <class Tup, class Func, std::size_t... Idx, class... Args>
			void __Proc(std::index_sequence<Idx...>, const Func& func, Args&&... args) {
				func(
					std::tuple_element_t<Idx,Tup>()...,
					std::forward<Args>(args)...
				);
			}
			template <class Tup, class Func, class... Args>
			void _Proc(const Func&, lubee::IConst<-1>, Args&&...) {}
			template <class Tup, class Func, int N, class... Args>
			void _Proc(const Func& func, lubee::IConst<N>, Args&&... args) {
				using Cur = std::tuple_element_t<N, Tup>;
				__Proc<Cur>(std::make_index_sequence<std::tuple_size<Cur>{}>(), func, std::forward<Args>(args)...);
				_Proc<Tup>(func, lubee::IConst<N-1>(), std::forward<Args>(args)...);
			}
			template <class Tup, class Func, class... Args>
			void Proc(const Func& func, Args&&... args) {
				_Proc<Tup>(
					func,
					lubee::IConst<std::tuple_size<Tup>{}-1>(),
					std::forward<Args>(args)...
				);
			}
			struct CheckConvert {
				template <int To, class V, class CB0, class CB1>
				static void _Proc(const V& src, const CB0& cb0, const CB1& cb1) {
					using value_t = typename V::value_t;
					Array<value_t, To> ar = {};
					using Cmp = lubee::Arithmetic<V::size, To>;
					for(int i=V::size ; i<To ; i++)
						ar[i] = 0;
					cb0(ar);
					for(int i=0 ; i<Cmp::less ; i++)
						ar[i] = src[i];

					using vec2_t = typename V::template type_cn<To>;
					const vec2_t v0(cb1(src)),
								v1(cb1(src.asInternal()));
					constexpr auto Th = lubee::Threshold<value_t>(0.1, 0);
					ASSERT_LE(AbsMax(ar - v0), Th);
					ASSERT_LE(AbsMax(ar - v1), Th);
				}
				template <int To, class V>
				void operator()(lubee::IConst<To>, const V& src) const {
					_Proc<To>(
						src,
						[](auto&) {},
						[](const auto& src) {
							return src.template convert<To>();
						}
					);
				}
				template <int To, int Pos, class V, class I, ENABLE_IF((Pos>=To))>
				void operator()(lubee::IConst<To>, lubee::IConst<Pos>, const V&, const I&) const {}
				template <int To, int Pos, class V, class I, ENABLE_IF((Pos<To))>
				void operator()(lubee::IConst<To>, lubee::IConst<Pos>, const V& src, const I& val) const {
					_Proc<To>(
						src,
						[val](auto& ar) {
							ar[Pos] = val;
						},
						[val](const auto& src){
							return src.template convertI<To,Pos>(val);
						}
					);
				}
			};
		}
		TYPED_TEST(Vector, ConvertI) {
			using ConvertSeq_t = lubee::seq::ExpandTypes_t2<
				std::tuple,
				std::tuple<
					lubee::seq::Range_t<2,5>,
					lubee::seq::Range_t<2,5>
				>
			>;
			USING(value_t);
			const value_t val = this->mt().template getUniform<value_t>();
			const auto v = this->makeRVec(DefaultRange<value_t>);
			ASSERT_NO_FATAL_FAILURE((Proc<ConvertSeq_t>(CheckConvert(), v, val)));
		}
		TYPED_TEST(Vector, Convert) {
			using ConvertSeq_t = lubee::seq::ExpandTypes_t2<
				std::tuple,
				std::tuple<
					lubee::seq::Range_t<2,5>
				>
			>;
			USING(value_t);
			const auto v = this->makeRVec(DefaultRange<value_t>);
			ASSERT_NO_FATAL_FAILURE((Proc<ConvertSeq_t>(CheckConvert(), v)));
		}
		namespace {
			template <class T, class V, ENABLE_IF(!(lubee::HasIndex_t<T,int>{}))>
			bool IsInRange(const T& val, const V& tMin, const V& tMax) {
				return val>= tMin && val<=tMax;
			}
			template <class T, class V, ENABLE_IF((lubee::HasIndex_t<T,int>{}))>
			bool IsInRange(const T& val, const V& tMin, const V& tMax) {
				for(auto& v : val) {
					if(!IsInRange(v, tMin, tMax))
						return false;
				}
				return true;
			}
		}
		TYPED_TEST(Vector, Saturation) {
			USING(value_t);
			USING(array_t);
			USING(vec_t);

			const auto vec = this->makeRVec(DefaultRange<value_t>);
			const auto range = lubee::random::GenRange<value_t>(this->mt().template getUniformF<value_t>());
			array_t ar(vec);
			for(auto& a : ar)
				a = Saturate(a, range.from, range.to);
			constexpr auto Th = lubee::Threshold<value_t>(0.2, 0);
			const vec_t v0 = vec.saturation(range.from, range.to),
						v1 = vec.asInternal().saturation(range.from, range.to);
			ASSERT_LE(AbsMax(ar-v0), Th);
			ASSERT_LE(AbsMax(ar-v1), Th);
			ASSERT_TRUE(IsInRange(v0, range.from, range.to));
			ASSERT_TRUE(IsInRange(v1, range.from, range.to));
		}
		namespace {
			struct Check_Equality {
				template <int Pos, class V>
				void operator()(lubee::IConst<Pos>, const V& src) const {
					Array<typename V::value_t, V::size> ar(src);
					for(int i=0 ; i<V::size ; i++)
						ar[i] = ar[Pos];
					auto tmp = src.asInternal();
					tmp.template makeEquality<Pos>();

					using value_t = typename V::value_t;
					constexpr auto Th = lubee::Threshold<value_t>(0.1, 0);
					ASSERT_LE(AbsMax(ar-V(tmp)), Th);
				}
			};
			struct Check_MaskL {
				template <int Pos, class V>
				void operator()(lubee::IConst<Pos>, const V& src) const {
					Array<typename V::value_t, V::size> ar(src);
					for(int i=0 ; i<Pos+1 ; i++)
						ar[i] = 0;
					auto tmp0 = src.asInternal();
					auto tmp = tmp0.template maskL<Pos>();

					using value_t = typename V::value_t;
					constexpr auto Th = lubee::Threshold<value_t>(0.1, 0);
					ASSERT_LE(AbsMax(ar-V(tmp)), Th);
				}
			};
			struct Check_MaskH {
				template <int Pos, class V>
				void operator()(lubee::IConst<Pos>, const V& src) const {
					Array<typename V::value_t, V::size> ar(src);
					for(int i=Pos+1 ; i<V::size ; i++)
						ar[i] = 0;
					auto tmp0 = src.asInternal();
					auto tmp = tmp0.template maskH<Pos>();

					using value_t = typename V::value_t;
					constexpr auto Th = lubee::Threshold<value_t>(0.1, 0);
					ASSERT_LE(AbsMax(ar-V(tmp)), Th);
				}
			};
			struct Check_PickAt {
				template <int Pos, class V>
				void operator()(lubee::IConst<Pos>, const V& src) const {
					using value_t = typename V::value_t;
					const Array<value_t, V::size> ar(src);
					constexpr auto Th = lubee::Threshold<value_t>(0.1, 0);
					ASSERT_LE(std::abs(ar[Pos] - src.asInternal().template pickAt<Pos>()), Th);
				}
			};
			struct Check_InitAt {
				template <int Pos, class V, class I>
				void operator()(lubee::IConst<Pos>, const V& src, const I& val) const {
					Array<typename V::value_t, V::size> ar;
					for(auto& a : ar)
						a = 0;
					ar[Pos] = val;
					auto tmp = src.asInternal();
					tmp.template initAt<Pos>(val);
				}
			};
		}
		TYPED_TEST(Vector, InternalFunc) {
			USING(vec_t);
			using Types = lubee::seq::ExpandTypes_t2<
				std::tuple,
				std::tuple<
					lubee::seq::Range_t<0, vec_t::size-1>
				>
			>;
			USING(value_t);
			const auto v = this->makeRVec(DefaultRange<value_t>);
			// 各種ルーチンのチェック
			ASSERT_NO_FATAL_FAILURE(Proc<Types>(Check_Equality(), v));
			ASSERT_NO_FATAL_FAILURE(Proc<Types>(Check_MaskL(), v));
			ASSERT_NO_FATAL_FAILURE(Proc<Types>(Check_MaskH(), v));
			ASSERT_NO_FATAL_FAILURE(Proc<Types>(Check_PickAt(), v));
			const value_t val = this->mt().template getUniform<value_t>();
			ASSERT_NO_FATAL_FAILURE(Proc<Types>(Check_InitAt(), v, val));
		}
		TYPED_TEST(Vector, Average) {
			USING(value_t);
			USING(vec_t);

			const auto vec = this->makeRVec(DefaultRange<value_t>);
			const value_t a0 = vec.average(),
							a1 = vec.asInternal().average();
			constexpr auto Th = lubee::Threshold<value_t>(0.7, 0);
			ASSERT_NEAR(a0, a1, Th);

			value_t avg = 0;
			for(const auto& a : vec)
				avg += a;
			avg /= vec_t::size;

			ASSERT_NEAR(avg, a0, Th);
		}
		TYPED_TEST(Vector, DotProduct) {
			USING(value_t);
			USING(array_t);
			USING(vec_t);

			// 配列で計算した場合と比較
			constexpr lubee::Range<value_t> range{1e2};
			const vec_t v0 = this->makeRVec(range),
						v1 = this->makeRVec(range);
			const array_t ar0(v0),
							ar1(v1);
			const value_t res0 = v0.dot(v1),
							res1 = v0.asInternal().dot(v1);
			const array_t ar2 = ar0 * ar1;
			value_t sum = 0;
			for(auto& a : ar2)
				sum += a;

			constexpr auto Th = lubee::Threshold<value_t>(0.8, 0);
			ASSERT_NEAR(sum, res0, Th);
			ASSERT_NEAR(sum, res1, Th);

			// ベクトルを反転したらDotProductの結果の符号も反転
			constexpr auto ThInv = lubee::Threshold<value_t>(0.1, 0);
			const vec_t vi = -v0;
			const value_t res2 = vi.dot(v1);
			ASSERT_NEAR(-res2, res0, ThInv);
		}
		TYPED_TEST(Vector, Absolute) {
			USING(value_t);
			USING(vec_t);
			constexpr auto range = DefaultRange<value_t>;
			const vec_t v = this->makeRVec(range);
			const vec_t va = v.absolute();
			for(int i=0 ; i<vec_t::size ; i++) {
				EXPECT_EQ(va[i], std::abs(v[i]));
			}
		}
		TYPED_TEST(Vector, GetMinMax) {
			USING(value_t);
			USING(vec_t);
			constexpr auto range = DefaultRange<value_t>;
			const vec_t v = this->makeRVec(range);
			value_t vMin=v[0],
					vMax=v[0];
			for(int i=1 ; i<vec_t::size ; i++) {
				vMin = std::min(vMin, v[i]);
				vMax = std::max(vMax, v[i]);
			}
			const value_t vMin0 = v.getMinValue(),
							vMax0 = v.getMaxValue();
			EXPECT_EQ(vMin0, vMin);
			EXPECT_EQ(vMax0, vMax);
		}
		TYPED_TEST(Vector, IsZero) {
			USING(value_t);
			USING(vec_t);
			constexpr auto Td = lubee::Threshold<value_t>(0.1, 1);
			const value_t Th = this->mt().template getUniform<value_t>({Td, 1e3});
			vec_t v;
			for(int i=0 ; i<vec_t::size ; i++)
				v[i] = this->mt().template getUniform<value_t>({-Th+Td, Th-Td});
			// 許容範囲値まではTrue
			EXPECT_TRUE(v.isZero(Th));
			// それ以上ならFalse
			const int idx = this->mt().template getUniform<int>({0,vec_t::size-1});
			v[idx] += Th*3;
			EXPECT_FALSE(v.isZero(Th));
		}
		TYPED_TEST(Vector, MinMax) {
			USING(value_t);
			USING(vec_t);
			USING(array_t);
			constexpr auto range = DefaultRange<value_t>;
			const vec_t v = this->makeRVec(range),
						vt = this->makeRVec(range);
			const vec_t vMin = v.getMin(vt),
						vMax = v.getMax(vt);
			array_t aMin,
					aMax;
			for(int i=0 ; i<vec_t::size ; i++) {
				aMin[i] = std::min(vt[i], v[i]);
				aMax[i] = std::max(vt[i], v[i]);
			}
			constexpr auto Th = lubee::Threshold<value_t>(0.1, 0);
			ASSERT_LE(AbsMax(aMin - vMin), Th);
			ASSERT_LE(AbsMax(aMax - vMax), Th);

			auto vMin0 = v,
				 vMax0 = v;
			vMin0.selectMin(vt);
			vMax0.selectMax(vt);
			ASSERT_LE(AbsMax(aMin - vMin0), Th);
			ASSERT_LE(AbsMax(aMax - vMax0), Th);
		}
		TYPED_TEST(Vector, SumUp) {
			USING(value_t);
			USING(array_t);
			constexpr auto range = DefaultRange<value_t>;
			const auto vec = this->makeRVec(range);
			const array_t ar(vec);

			value_t sum0 = 0;
			for(int i=0 ; i<array_t::size ; i++)
				sum0 += ar[i];
			const value_t sum1 = (vec+0).sumUp();
			constexpr auto Th = lubee::Threshold<value_t>(0.8, 0);
			ASSERT_LE(sum0-sum1, Th);
		}
		// 内部表現による演算結果の差異をチェック
		TYPED_TEST(Vector, Register) {
			USING(value_t);
			auto& mt = this->mt();
			constexpr auto range = DefaultRange<value_t>;
			auto vec = this->makeRVec(range);
			constexpr int size = decltype(vec)::size;
			typename TestFixture::array_t	raw;
			for(int i=0 ; i<size ; i++)
				raw.m[i] = vec.m[i];
			const auto chk = [&](auto dist){
				for(int i=0 ; i<size ; i++) {
					if(TestFixture::integral)
						ASSERT_NEAR(raw.m[i], vec.m[i], dist);
					else
						ASSERT_FLOAT_EQ(raw.m[i], vec.m[i]);
				}
			};
			const auto rdi = mt.template getUniformF<int>();
			const auto t = this->makeRVecNZ(value_t(TestFixture::integral ? 1 : 1e-4), {-1e4, 1e4});
			// ランダムで四則演算
			switch(rdi({0,3})) {
				case 0:
					vec += t;
					raw += t;
					// 結果が一致することを確認
					ASSERT_NO_FATAL_FAILURE(chk(0));
					break;
				case 1:
					vec -= t;
					raw -= t;
					ASSERT_NO_FATAL_FAILURE(chk(0));
					break;
				case 2:
					vec *= t;
					raw *= t;
					ASSERT_NO_FATAL_FAILURE(chk(0));
					break;
				case 3:
					vec /= t;
					raw /= t;
					ASSERT_NO_FATAL_FAILURE(chk(0));
					break;
				default:
					__assume(false)
			}
		}
		// 要素比較チェック
		TYPED_TEST(Vector, Compare) {
			USING(value_t);
			constexpr auto range = DefaultRange<value_t>;
			const auto v0 = this->makeRVec(range),
						v1 = this->makeRVec(range);
			// 自身と比較したら==はtrueになる
			ASSERT_EQ(v0, v0);
			// ==と!=は正反対の結果になる
			ASSERT_NE(v0==v1, v0!=v1);

			USING(vec_t);
			const auto randNumNZ = [this]{
				auto ret = this->mt().template getUniform<value_t>();
				if(std::abs(ret) < 1)
					ret += 10;
				return ret;
			};
			// 自身を足し合わせれば値が変わる
			vec_t tv0 = v0 + v0;
			ASSERT_NE(tv0, v0);
			for(int i=0 ; i<vec_t::size ; i++) {
				tv0 = v0;
				ASSERT_EQ(tv0, v0);
				// 要素を1つだけいじっても==がfalseになる
				tv0.m[i] += randNumNZ();
				ASSERT_NE(tv0, v0);
			}
		}
		// 論理演算チェック
		TYPED_TEST(Vector, Logical) {
			USING(value_t);
			USING(vec_t);
			constexpr auto range = DefaultRange<value_t>;
			const auto v0 = this->makeRVec(range);
			ASSERT_EQ(v0.asInternal(), v0);
			// 自身とのAndは元と同じ
			ASSERT_EQ(v0&v0, v0);
			// 自身とのOrも同じ
			ASSERT_EQ(v0|v0, v0);
			// 自身とのXorは必ずゼロになる
			ASSERT_EQ(vec_t(0), v0^v0);

			// (整数値のみ)
			if(TestFixture::integral) {
				const vec_t one = vec_t::I::One();
				// ~0とのAndは元と同じ
				ASSERT_EQ(v0, v0&one);
				// ~0とのOrは~0になる
				ASSERT_EQ(one, v0|one);
				// ~0のXorで2回反転させれば元と同じになる
				ASSERT_EQ(v0, (v0^one^one));
			}
		}
		// 四則演算のチェック
		TYPED_TEST(Vector, Arithmetic) {
			USING(value_t);
			constexpr auto range = DefaultRange<value_t>;
			const auto v0 = this->makeRVec(range),
						v1 = this->makeRVec(range);
			const auto check = [](const auto& v0, const auto& v1, const auto& op, const auto& ope) {
				using v_t = std::decay_t<decltype(v0)>;
				{
					// 代入後の値チェック
					v_t tv0 = v0,
						tv1 = v1;
					ASSERT_EQ(v0, tv0);
					ASSERT_EQ(v1, tv1);

					// 三項演算をどこにも代入しなければ値は変わらない
					op(v0, v1);
					ASSERT_EQ(tv0, v0);
					ASSERT_EQ(tv1, v1);
				}
				{
					// 二項演算と三項演算は一緒
					v_t tv0 = op(v0, v1),
						tv1 = v0;
					ope(tv1, v1);
					ASSERT_EQ(tv0, tv1);
					// 複数の三項演算を二項演算に分けても結果は一緒
					tv0 = op(op(v0, v0), v1);
					tv1 = v0;
					ope(tv1, v0);
					ope(tv1, v1);
					ASSERT_EQ(tv0, tv1);
				}
			};
			#define DEF_TEST(op) \
				ASSERT_NO_FATAL_FAILURE( \
					check( \
						v0, v1, \
						[](const auto& p0, const auto& p1){ return p0 op p1; }, \
						[](auto& p0, const auto& p1){ return p0 op##= p1; } \
					) \
				);
			DEF_TEST(+)
			DEF_TEST(-)
			DEF_TEST(*)
			// Divideの場合は0除算を避ける
			constexpr auto Th = lubee::Threshold<value_t>(0.4, 1);
			if(!HasZero(v0.m, Th) && !HasZero(v1.m, Th)) {
				DEF_TEST(/)
			}
			#undef DEF_TEST
		}
		TYPED_TEST(Vector, Serialization) {
			lubee::CheckSerialization(this->makeRVec());
		}
	}
}
