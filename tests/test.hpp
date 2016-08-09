#pragma once
#include "../random.hpp"
#include "../matrix.hpp"
#include "../random/vector.hpp"
#include "../random/matrix.hpp"
#include "../random/range.hpp"
#include "../random/quaternion.hpp"
#include "../ulps.hpp"
#include "../meta/check_macro.hpp"
#include <gtest/gtest.h>

namespace frea {
	namespace test {
		template <class T, class TV>
		bool HasZero(const T& t, const TV& th) {
			for(auto& tm : t) {
				if(std::abs(tm) < th)
					return true;
			}
			return false;
		}
		template <class T>
		auto DiffSum(const T& t0, const T& t1) {
			static_assert(countof(t0) == countof(t1), "");
			auto sum = std::abs(t0[0] - t1[0]);
			for(int i=1 ; i<int(countof(t0)) ; i++)
				sum += std::abs(t0[i] - t1[i]);
			return sum;
		}
		template <class T, class RDF>
		void _Fill(T& dst, const RDF& rdf, ...) {
			dst = rdf();
		}
		template <class T, class RDF>
		void _Fill(T& dst, const RDF& rdf, std::decay_t<decltype(std::declval<T>()[0])>*) {
			for(auto& d : dst.m)
				_Fill(d, rdf, nullptr);
		}
		template <class T, class RDF>
		void Fill(T& dst, const RDF& rdf) {
			_Fill(dst, rdf, nullptr);
		}

		template <class T, class OP>
		auto _Op(const T& t, const OP&, ...) {
			return t;
		}
		template <class T, class OP>
		auto _Op(const T& t, const OP& op, std::decay_t<decltype(std::declval<T>()[0])>*) {
			auto ret = _Op(t.m[0], op, nullptr);
			for(int i=0 ; i<int(countof(t.m)) ; i++)
				ret = op(ret, _Op(t.m[i], op, nullptr));
			return ret;
		}
		template <class T, class OP>
		auto Op(const T& t, const OP& op) {
			return _Op(t, op, nullptr);
		}

		template <class T>
		auto SumUp(const T& t) {
			return Op(t, [](const auto& t0, const auto& t1){ return t0+t1; });
		}
		template <class T>
		auto AbsMax(const T& t) {
			return Op(t, [](const auto& t0, const auto& t1){ return std::max(t0, std::abs(t1)); });
		}

		template <class T>
		auto LowerSize(std::decay_t<decltype(T::size)>*) -> IConst<T::size>;
		template <class T>
		auto LowerSize(...) -> IConst<1>;

		template <class T, int N>
		struct Array {
			constexpr static int size = N,
								lower_size = decltype(LowerSize<T>(nullptr))::value;
			T	m[size];

			#define DEF_OP(op) \
				template <class V> \
				Array& operator op##= (const V& r) { \
					return *this = *this op r; \
				} \
				template <class V, \
						ENABLE_IF((HasIndex_t<V,int>{}))> \
				Array operator op (const V& r) const { \
					static_assert(V::size==size, ""); \
					Array ret; \
					for(int i=0 ; i<size ; i++) \
						ret[i] = m[i] op r[i]; \
					return ret; \
				} \
				template <class V, \
						ENABLE_IF(!(HasIndex_t<V,int>{}))> \
				Array operator op (const V& r) const { \
					Array ret; \
					for(int i=0 ; i<size ; i++) \
						ret[i] = m[i] op r; \
					return ret; \
				}
			DEF_OP(+)
			DEF_OP(-)
			DEF_OP(/)
			DEF_OP(*)
			#undef DEF_OP
			Array() = default;
			template <class T2>
			Array(const Array<T2,size>& a): m(a.m) {}
			template <class V,
					 ENABLE_IF((HasIndex_t<V,int>{}))>
			Array(const V& v) {
				*this = v;
			}
			template <class V,
					 ENABLE_IF((HasIndex_t<V,int>{}))>
			Array& operator = (const V& r) {
				for(int i=0 ; i<size ; i++)
					m[i] = r[i];
				return *this;
			}

			T& operator [](const int n) noexcept {
				return m[n];
			}
			const T& operator [](const int n) const noexcept {
				return m[n];
			}
			template <class V>
			bool operator < (const V& v) const {
				static_assert(size==V::size, "");
				return std::lexicographical_compare(
							m, m+size,
							v.m, v.m+size,
							std::less<>()
						);
			}
			template <class V>
			bool near(const T& th, const V& v) const {
				static_assert(size==V::size, "");
				for(int i=0 ; i<size ; i++) {
					if(std::abs(m[i] - v.m[i]) > th)
						return false;
				}
				return true;
			}
			template <class V>
			bool operator == (const V& v) const {
				return near(T(0), v);
			}
			template <class V>
			bool operator != (const V& v) const {
				return !(operator ==(v));
			}
		};
		class Random : public ::testing::Test {
			private:
				RandomMT	_mt;
			public:
				Random(): _mt(RandomMT::Make<4>()) {}
				auto& mt() {
					return _mt;
				}
		};
		/*!
			\tparam tuple<レジスタタイプ, 要素数, アラインメント>
		*/
		template <class T>
		class RVector : public Random {
			public:
				using reg_t = std::tuple_element_t<0,T>;
				constexpr static int size = std::tuple_element_t<1,T>::value;
				constexpr static bool align = std::tuple_element_t<2,T>::value;
				using vec_t = SVec_t<reg_t, size, align>;
				using value_t = typename vec_t::value_t;
				constexpr static bool integral = std::is_integral<value_t>::value;
				using array_t = Array<value_t, size>;
			private:
				using RFunc = std::function<value_t (const Range<value_t>&)>;
				RFunc	_rdf;
			public:
				RVector():
					_rdf(this->mt().template getUniformF<value_t>())
				{}
				auto& rdf() noexcept { return _rdf; }
				auto makeRVec(const Range<value_t>& r=random::DefaultVecRange<value_t>) {
					return random::GenVecRange<vec_t>(_rdf, r);
				}
				auto makeRVecNZ(const value_t& th, const Range<value_t>& r=random::DefaultVecRange<value_t>) {
					return random::GenVecCnd<vec_t>(_rdf,
							[th](const auto& v){
								return !HasZero(v.m, th);
							}, r);
				}
		};
		template <class T, int M, int N>
		struct ArrayM : Array<Array<T,N>, M> {
			using base_t = Array<Array<T,N>, M>;
			using base_t::base_t;
			template <class V2>
			auto mul(const V2& v, std::false_type) const {
				return base_t::operator * (v);
			}
			template <class A2>
			auto mul(const A2& m, std::true_type) const {
				static_assert(N==A2::lower_size, "invalid operation");
				ArrayM<T, M, A2::lower_size> ret;
				for(int i=0 ; i<M ; i++) {
					for(int j=0 ; j<A2::lower_size ; j++) {
						auto& dst = ret[i][j];
						dst = 0;
						for(int k=0 ; k<base_t::size ; k++) {
							dst += (*this)[i][k] * m[k][j];
						}
					}
				}
				return ret;
			}
			template <class A2>
			auto operator * (const A2& m) const {
				return mul(m, HasIndex_t<A2,int>());
			}
			template <class A2>
			ArrayM& operator *= (const A2& m) {
				return *this = *this * m;
			}
		};
		/*!
			\tparam tuple<レジスタタイプ, 要素数M, 要素数N, アラインメント>
		*/
		template <class T>
		class RMatrix : public RVector<seq::StripAt_t<T,1>> {
			public:
				using base_t = RVector<seq::StripAt_t<T,1>>;
				constexpr static int dim_m = std::tuple_element_t<1,T>::value,
									dim_n = std::tuple_element_t<2,T>::value;
				using reg_t = typename base_t::reg_t;
				using value_t = typename base_t::value_t;
				constexpr static auto align = base_t::align;
				using mat_t = SMat_t<reg_t, dim_m, dim_n, align>;
				using array_t = ArrayM<value_t, dim_m, dim_n>;
			public:
				auto makeRMat(const Range<value_t>& r=random::DefaultVecRange<value_t>) {
					return random::GenMatRange<mat_t>(base_t::rdf(), r);
				}
		};

		template <class T>
		constexpr T ThresholdULPs;
		template <>
		constexpr auto ThresholdULPs<float> = ulps::Diff_C<float>(0.f, 3e-4f);
		template <>
		constexpr auto ThresholdULPs<double> = ulps::Diff_C<double>(0.f, 3e-4);

		template <class T>
		constexpr T RangeV;
		template <>
		constexpr auto RangeV<float> = 1e4f;
		template <>
		constexpr auto RangeV<double> = 1e16;

		template <class... Ts0, class... Ts1>
		auto ConcatTypes(::testing::Types<Ts0...>, ::testing::Types<Ts1...>) -> ::testing::Types<Ts0..., Ts1...>;
		template <class T0, class T1>
		using ConcatTypes_t = decltype(ConcatTypes(std::declval<T0>(), std::declval<T1>()));

		template <class... Ts>
		auto ToTestTypes(std::tuple<Ts...>) -> ::testing::Types<Ts...>;
		template <class T>
		using ToTestTypes_t = decltype(ToTestTypes(std::declval<T>()));

		// 整数ベクトル
		using ITypes_t = seq::ExpandTypes_t2<
				std::tuple,
				std::tuple<
					std::tuple<__m128i>,
					seq::Range_t<2,5>,
					seq::BoolSeq_t
				>
		>;
		// 浮動小数点数ベクトル
		using FTypes_t = seq::ExpandTypes_t2<
				std::tuple,
				std::tuple<
					std::tuple<__m128, __m128d>,
					seq::Range_t<2,5>,
					seq::BoolSeq_t
				>
			>;
		// 各要素数固有の関数テスト用
		template <int N>
		using TypesD_t = seq::ExpandTypes_t2<
				std::tuple,
				std::tuple<
					std::tuple<__m128, __m128d>,
					std::tuple<IConst<N>>,
					std::tuple<BConst<false>>
				>
		>;
	}
}
