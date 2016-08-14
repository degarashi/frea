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
		//! デバッガでチェックしやすいように要素を等差数列で初期化
		template <class T, std::size_t... Idx>
		T ArithmeticSequence(std::index_sequence<Idx...>) {
			return T(static_cast<typename T::value_t>(Idx)...);
		}

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

		template <class T, class OP, ENABLE_IF(!(HasIndex_t<T,int>{}))>
		auto Op(const T& t, const OP&) {
			return t;
		}
		template <class T, class OP, ENABLE_IF((HasIndex_t<T,int>{}))>
		auto Op(const T& t, const OP& op) {
			auto ret = Op(t.m[0], op);
			for(int i=0 ; i<int(countof(t.m)) ; i++)
				ret = op(ret, Op(t.m[i], op));
			return ret;
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
		class Array {
			public:
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
				using elem_t = std::tuple_element_t<0,T>;
				constexpr static int size = std::tuple_element_t<1,T>::value;
				constexpr static bool align = std::tuple_element_t<2,T>::value;
				using vec_t = Vec_t<elem_t, size, align>;
				using value_t = typename vec_t::value_t;
				constexpr static bool integral = std::is_integral<value_t>::value;
				using array_t = Array<value_t, size>;
			public:
				auto makeRVec() {
					return random::GenVec<vec_t>(this->mt().template getUniformF<value_t>());
				}
				auto makeRVec(const Range<value_t>& r) {
					return random::GenVec<vec_t>(this->mt().template getUniformF<value_t>(r));
				}
				auto makeRVecNZ(const value_t& th, const Range<value_t>& r) {
					auto rd = this->mt().template getUniformF<value_t>(r);
					return random::GenVecCnd<vec_t>(rd,
							[th](const auto& v){
								return !HasZero(v.m, th);
							});
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
			void transpose() {
				for(int i=0 ; i<base_t::size ; i++) {
					for(int j=i+1 ; j<base_t::lower_size ; j++) {
						std::swap((*this)[i][j], (*this)[j][i]);
					}
				}
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
				using value_t = typename base_t::value_t;
				constexpr static auto align = base_t::align;
				using elem_t = typename base_t::elem_t;
				using mat_t = Mat_t<elem_t, dim_m, dim_n, align>;
				using array_t = ArrayM<value_t, dim_m, dim_n>;
			public:
				auto makeRMat(const Range<value_t>& r) {
					auto rd = this->mt().template getUniformF<value_t>(r);
					return random::GenMat<mat_t>(rd);
				}
		};

		template <class T, ENABLE_IF(std::is_floating_point<T>{})>
		constexpr T Threshold(const T& fr, const T&) noexcept {
			return std::numeric_limits<T>::epsilon()*fr;
		}
		template <class T, ENABLE_IF(std::is_integral<T>{})>
		constexpr T Threshold(const T&, const T& ia) noexcept {
			return ia;
		}
		template <class T>
		constexpr bool IsZero(const T& val, const T& th) noexcept {
			return std::abs(val) < th;
		}

		template <class... Ts0, class... Ts1>
		auto ConcatTypes(::testing::Types<Ts0...>, ::testing::Types<Ts1...>) -> ::testing::Types<Ts0..., Ts1...>;
		template <class T0, class T1>
		using ConcatTypes_t = decltype(ConcatTypes(std::declval<T0>(), std::declval<T1>()));

		template <class... Ts>
		auto ToTestTypes(std::tuple<Ts...>) -> ::testing::Types<Ts...>;
		template <class T>
		using ToTestTypes_t = decltype(ToTestTypes(std::declval<T>()));

		namespace types {
			using FElem_t = std::tuple<float, double>;
			using IElem_t = std::tuple<int32_t>;
			using Elem_t = seq::TupleCat_t<FElem_t, IElem_t>;

			#if SSE>=2
				using FReg_t = std::tuple<__m128, __m128d>;
				using IReg_t = std::tuple<__m128i>;
			#else
				using FReg_t = std::tuple<>;
				using IReg_t = std::tuple<>;
			#endif
			using Reg_t = seq::TupleCat_t<FReg_t, IReg_t>;

			using Float_t = seq::TupleCat_t<Elem_t, Reg_t>;
			using Int_t = seq::TupleCat_t<IElem_t, IReg_t>;
			using Value_t = seq::TupleCat_t<Float_t, Int_t>;

			template <class E>
			using VectorRange_t = seq::ExpandTypes_t2<
				std::tuple,
				std::tuple<
					E,
					seq::Range_t<2,5>,
					seq::BoolSeq_t
				>
			>;

			// 浮動小数点数ベクトル
			// 各要素数固有の関数テスト用
			template <class E, int N>
			using VectorRangeD_t = seq::ExpandTypes_t2<
				std::tuple,
				std::tuple<
					E,
					std::tuple<IConst<N>>,
					std::tuple<BConst<false>>
				>
			>;

			template <class R, class S, class A>
			using SquareMat_t = std::tuple<R, S, S, A>;
			template <class E, int M, int N>
			using MatrixRange_t = seq::ExpandTypes_t2<
				std::tuple,
				std::tuple<
					E,
					seq::Range_t<2,M>,
					seq::Range_t<2,N>,
					seq::BoolSeq_t
				>
			>;

			using QTypes_t = seq::ExpandTypes_t2<
				std::tuple,
				std::tuple<
					FElem_t,
					seq::BoolSeq_t
				>
			>;
			using QTypes = ToTestTypes_t<QTypes_t>;
		}
	}
}
