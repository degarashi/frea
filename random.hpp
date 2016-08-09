#pragma once
#include "meta/enable_if.hpp"
#include "range.hpp"
#include <algorithm>
#include <array>
#include <functional>

namespace frea {
	template <class MT>
	class Random {
		private:
			using MT_t = MT;
			MT_t	_mt;

			template <class T>
			using Lim = std::numeric_limits<T>;
			template <class T>
			struct Dist_Int {
				using uniform_t = std::uniform_int_distribution<T>;
				constexpr static Range<T> DefaultRange{Lim<T>::lowest(), Lim<T>::max()},
												NumericRange = DefaultRange;
			};
			template <class T>
			struct Dist_Float  {
				using uniform_t = std::uniform_real_distribution<T>;
				constexpr static Range<T> DefaultRange{T(0), T(1)},
												NumericRange{Lim<T>::lowest()/2, Lim<T>::max()/2};
			};
			template <class T, ENABLE_IF((std::is_integral<T>{}))>
			static Dist_Int<T> DetectDist();
			template <class T, ENABLE_IF((std::is_floating_point<T>{}))>
			static Dist_Float<T> DetectDist();
			template <class T>
			using Dist_t = decltype(DetectDist<T>());
		public:
			template <class T>
			class RObj {
				private:
					Random&			_rd;
					const Range<T>	_range;
				public:
					RObj(Random& rd, const Range<T>& r):
						_rd(rd),
						_range(r)
					{}
					RObj(const RObj&) = default;
					//! 実行時範囲指定
					auto operator ()(const Range<T>& r) const {
						return _rd.template getUniform<T>(r);
					}
					//! デフォルト範囲指定
					auto operator ()() const {
						return (*this)(_range);
					}
			};
		public:
			Random(MT_t mt) noexcept: _mt(std::move(mt)) {}
			Random(const Random&) = delete;
			Random(Random&& t) noexcept {
				*this = std::move(t);
			}
			Random& operator = (Random&& t) noexcept {
				swap(t);
				return *this;
			}
			void swap(Random& t) noexcept {
				std::swap(_mt, t._mt);
			}
			MT_t& refMt() noexcept {
				return _mt;
			}
			//! 一様分布
			/*! floating-pointの場合は0から1の乱数
				integerの場合は範囲指定なしnumeric_limits<T>::min() -> max()の乱数 */
			template <class T>
			T getUniform() {
				constexpr auto R = Dist_t<T>::DefaultRange;
				return getUniform<T>(R);
			}
			//! 指定範囲の一様分布(in range)
			template <class T>
			T getUniform(const Range<T>& range) {
				return typename Dist_t<T>::uniform_t(range.from, range.to)(_mt);
			}
			//! 一様分布を返すファンクタを作成
			template <class T>
			auto getUniformF(const Range<T>& r) noexcept {
				return RObj<T>(*this, r);
			}
			//! 一様分布を返すファンクタを作成
			template <class T>
			auto getUniformF() noexcept {
				constexpr auto R = Dist_t<T>::DefaultRange;
				return RObj<T>(*this, R);
			}
			//! 指定範囲の一様分布(vmax)
			template <class T>
			T getUniformMax(const T& vmax) {
				return getUniform<T>({Lim<T>::lowest(), vmax});
			}
			//! 指定範囲の一様分布(vmin)
			template <class T>
			T getUniformMin(const T& vmin) {
				return getUniform<T>({vmin, Lim<T>::max()});
			}

			template <std::size_t NumSeed>
			static Random<MT_t> Make() {
				std::random_device rnd;
				std::array<std::seed_seq::result_type, NumSeed> seed;
				std::generate(seed.begin(), seed.end(), std::ref(rnd));
				std::seed_seq seq(seed.begin(), seed.end());
				return std::move(Random<MT_t>(MT_t{seq}));
			}
	};
	using RandomMT = Random<std::mt19937>;
	using RandomMT64 = Random<std::mt19937_64>;
}
