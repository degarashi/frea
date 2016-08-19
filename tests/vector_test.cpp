#include "test.hpp"
#include "../compiler_macro.hpp"

namespace frea {
	namespace test {
		// Int, Float両方のテストケース
		template <class T>
		using Vector = RVector<T>;
		using Types = ToTestTypes_t<types::VectorRange_t<types::Reg_t>>;
		TYPED_TEST_CASE(Vector, Types);

		namespace {
			template <class T>
			constexpr Range<T> DefaultRange{-1e3, 1e3};
		}
		TYPED_TEST(Vector, Minus) {
			using value_t = typename TestFixture::value_t;
			using array_t = typename TestFixture::array_t;
			using vec_t = typename TestFixture::vec_t;

			// operator -
			constexpr auto range = DefaultRange<value_t>;
			auto vec = this->makeRVec(range);
			// 地道に要素の符号反転した結果と比べる
			array_t ar(vec);
			for(auto& a : ar.m)
				a = -a;
			vec = -vec;
			constexpr auto Th = Threshold<value_t>(0.1, 0);
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
			void _Proc(const Func&, IConst<-1>, Args&&...) {}
			template <class Tup, class Func, int N, class... Args>
			void _Proc(const Func& func, IConst<N>, Args&&... args) {
				using Cur = std::tuple_element_t<N, Tup>;
				__Proc<Cur>(std::make_index_sequence<std::tuple_size<Cur>{}>(), func, std::forward<Args>(args)...);
				_Proc<Tup>(func, IConst<N-1>(), std::forward<Args>(args)...);
			}
			template <class Tup, class Func, class... Args>
			void Proc(const Func& func, Args&&... args) {
				_Proc<Tup>(
					func,
					IConst<std::tuple_size<Tup>{}-1>(),
					std::forward<Args>(args)...
				);
			}
			struct CheckConvert {
				template <int To, class V, class CB0, class CB1>
				static void _Proc(const V& src, const CB0& cb0, const CB1& cb1) {
					using value_t = typename V::value_t;
					Array<value_t, To> ar = {};
					cb0(ar);
					using Cmp = Arithmetic<V::size, To>;
					for(int i=0 ; i<Cmp::less ; i++)
						ar[i] = src[i];
					for(int i=V::size ; i<To ; i++)
						ar[i] = 0;

					using vec2_t = typename V::template type_cn<To>;
					const vec2_t v0(cb1(src)),
								v1(cb1(src.asInternal()));
					constexpr auto Th = Threshold<value_t>(0.1, 0);
					ASSERT_LE(AbsMax(ar - v0), Th);
					ASSERT_LE(AbsMax(ar - v1), Th);
				}
				template <int To, class V>
				void operator()(IConst<To>, const V& src) const {
					_Proc<To>(
						src,
						[](auto&) {},
						[](const auto& src) {
							return src.template convert<To>();
						}
					);
				}
				template <int To, int Pos, class V, class I, ENABLE_IF((Pos>=To))>
				void operator()(IConst<To>, IConst<Pos>, const V&, const I&) const {}
				template <int To, int Pos, class V, class I, ENABLE_IF((Pos<To))>
				void operator()(IConst<To>, IConst<Pos>, const V& src, const I& val) const {
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
			using ConvertSeq_t = seq::ExpandTypes_t2<
				std::tuple,
				std::tuple<
					seq::Range_t<2,5>,
					seq::Range_t<2,5>
				>
			>;
			using value_t = typename TestFixture::value_t;
			const value_t val = this->mt().template getUniform<value_t>();
			const auto v = this->makeRVec(DefaultRange<value_t>);
			ASSERT_NO_FATAL_FAILURE((Proc<ConvertSeq_t>(CheckConvert(), v, val)));
		}
		TYPED_TEST(Vector, Convert) {
			using ConvertSeq_t = seq::ExpandTypes_t2<
				std::tuple,
				std::tuple<
					seq::Range_t<2,5>
				>
			>;
			using value_t = typename TestFixture::value_t;
			const auto v = this->makeRVec(DefaultRange<value_t>);
			ASSERT_NO_FATAL_FAILURE((Proc<ConvertSeq_t>(CheckConvert(), v)));
		}
		namespace {
			template <class T, class V, ENABLE_IF(!(HasIndex_t<T,int>{}))>
			bool IsInRange(const T& val, const V& tMin, const V& tMax) {
				return val>= tMin && val<=tMax;
			}
			template <class T, class V, ENABLE_IF((HasIndex_t<T,int>{}))>
			bool IsInRange(const T& val, const V& tMin, const V& tMax) {
				for(auto& v : val) {
					if(!IsInRange(v, tMin, tMax))
						return false;
				}
				return true;
			}
		}
		TYPED_TEST(Vector, Saturation) {
			using value_t = typename TestFixture::value_t;
			using array_t = typename TestFixture::array_t;
			using vec_t = typename TestFixture::vec_t;

			const auto vec = this->makeRVec(DefaultRange<value_t>);
			const auto range = random::GenRange<value_t>(this->mt().template getUniformF<value_t>());
			array_t ar(vec);
			for(auto& a : ar)
				a = Saturate(a, range.from, range.to);
			constexpr auto Th = Threshold<value_t>(0.2, 0);
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
				void operator()(IConst<Pos>, const V& src) const {
					Array<typename V::value_t, V::size> ar(src);
					for(int i=0 ; i<V::size ; i++)
						ar[i] = ar[Pos];
					auto tmp = src.asInternal();
					tmp.template makeEquality<Pos>();

					using value_t = typename V::value_t;
					constexpr auto Th = Threshold<value_t>(0.1, 0);
					ASSERT_LE(AbsMax(ar-V(tmp)), Th);
				}
			};
			struct Check_MaskL {
				template <int Pos, class V>
				void operator()(IConst<Pos>, const V& src) const {
					Array<typename V::value_t, V::size> ar(src);
					for(int i=0 ; i<Pos+1 ; i++)
						ar[i] = 0;
					auto tmp0 = src.asInternal();
					auto tmp = tmp0.template maskL<Pos>();

					using value_t = typename V::value_t;
					constexpr auto Th = Threshold<value_t>(0.1, 0);
					ASSERT_LE(AbsMax(ar-V(tmp)), Th);
				}
			};
			struct Check_MaskH {
				template <int Pos, class V>
				void operator()(IConst<Pos>, const V& src) const {
					Array<typename V::value_t, V::size> ar(src);
					for(int i=Pos+1 ; i<V::size ; i++)
						ar[i] = 0;
					auto tmp0 = src.asInternal();
					auto tmp = tmp0.template maskH<Pos>();

					using value_t = typename V::value_t;
					constexpr auto Th = Threshold<value_t>(0.1, 0);
					ASSERT_LE(AbsMax(ar-V(tmp)), Th);
				}
			};
			struct Check_PickAt {
				template <int Pos, class V>
				void operator()(IConst<Pos>, const V& src) const {
					using value_t = typename V::value_t;
					const Array<value_t, V::size> ar(src);
					constexpr auto Th = Threshold<value_t>(0.1, 0);
					ASSERT_LE(std::abs(ar[Pos] - src.asInternal().template pickAt<Pos>()), Th);
				}
			};
			struct Check_InitAt {
				template <int Pos, class V, class I>
				void operator()(IConst<Pos>, const V& src, const I& val) const {
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
			using vec_t = typename TestFixture::vec_t;
			using Types = seq::ExpandTypes_t2<
				std::tuple,
				std::tuple<
					seq::Range_t<0, vec_t::size-1>
				>
			>;
			using value_t = typename TestFixture::value_t;
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
			using value_t = typename TestFixture::value_t;
			using vec_t = typename TestFixture::vec_t;

			const auto vec = this->makeRVec(DefaultRange<value_t>);
			const value_t a0 = vec.average(),
							a1 = vec.asInternal().average();
			constexpr auto Th = Threshold<value_t>(0.3, 0);
			ASSERT_NEAR(a0, a1, Th);

			value_t avg = 0;
			for(const auto& a : vec)
				avg += a;
			avg /= vec_t::size;

			ASSERT_NEAR(avg, a0, Th);
		}
		TYPED_TEST(Vector, DotProduct) {
			using value_t = typename TestFixture::value_t;
			using array_t = typename TestFixture::array_t;
			using vec_t = typename TestFixture::vec_t;

			// 配列で計算した場合と比較
			constexpr auto range = DefaultRange<value_t>;
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

			constexpr auto Th = Threshold<value_t>(0.6, 0);
			ASSERT_NEAR(sum, res0, Th);
			ASSERT_NEAR(sum, res1, Th);

			// ベクトルを反転したらDotProductの結果の符号も反転
			constexpr auto ThInv = Threshold<value_t>(0.1, 0);
			const vec_t vi = -v0;
			const value_t res2 = vi.dot(v1);
			ASSERT_NEAR(-res2, res0, ThInv);
		}
		TYPED_TEST(Vector, MinMax) {
			using value_t = typename TestFixture::value_t;
			using vec_t = typename TestFixture::vec_t;
			using array_t = typename TestFixture::array_t;
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
			constexpr auto Th = Threshold<value_t>(0.1, 0);
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
			using value_t = typename TestFixture::value_t;
			using array_t = typename TestFixture::array_t;
			constexpr auto range = DefaultRange<value_t>;
			const auto vec = this->makeRVec(range);
			const array_t ar(vec);

			value_t sum0 = 0;
			for(int i=0 ; i<array_t::size ; i++)
				sum0 += ar[i];
			const value_t sum1 = (vec+0).sumUp();
			constexpr auto Th = Threshold<value_t>(0.4, 0);
			ASSERT_LE(sum0-sum1, Th);
		}
		// 内部表現による演算結果の差異をチェック
		TYPED_TEST(Vector, Register) {
			using value_t = typename TestFixture::value_t;
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
			using value_t = typename TestFixture::value_t;
			constexpr auto range = DefaultRange<value_t>;
			const auto v0 = this->makeRVec(range),
						v1 = this->makeRVec(range);
			// 自身と比較したら==はtrueになる
			ASSERT_EQ(v0, v0);
			// ==と!=は正反対の結果になる
			ASSERT_NE(v0==v1, v0!=v1);

			using vec_t = typename TestFixture::vec_t;
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
			using value_t = typename TestFixture::value_t;
			using vec_t = typename TestFixture::vec_t;
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
			using value_t = typename TestFixture::value_t;
			constexpr auto range = DefaultRange<value_t>;
			const auto v0 = this->makeRVec(range),
						v1 = this->makeRVec(range);
			const auto check = [](auto v0, auto v1, const auto& op, const auto& ope) {
				// 代入後の値チェック
				auto tv0 = v0;
				const auto tv1 = v1;
				ASSERT_EQ(v0, tv0);
				ASSERT_EQ(v1, tv1);

				// 三項演算をどこにも代入しなければ値は変わらない
				op(v0, v1);
				ASSERT_EQ(tv0, v0);
				ASSERT_EQ(tv1, v1);
				// 二項演算と三項演算は一緒
				tv0 = op(v0, v1);
				ope(v0, v1);
				ASSERT_EQ(tv0, v0);
				// 複数の三項演算を二項演算に分けても結果は一緒
				v0 = tv0;
				v0 = op(op(v0, v0), v1);
				ope(tv0, tv0);
				ope(tv0, v1);
				ASSERT_EQ(tv0, v0);
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
			if(!HasZero(v1.m, Threshold<value_t>(0.4, 1))) {
				DEF_TEST(/)
			}
			#undef DEF_TEST
		}
	}
}
