#ifndef _SINGULAR_SVD_H
#define _SINGULAR_SVD_H

#include "singular/Matrix.h"
#include "singular/Reflector.h"
#include "singular/Rotator.h"

#include <iostream>
#include <tuple>

namespace singular {

	/**
	 * Namespace for singular value decomposition.
	 *
	 * @tparam M
	 *     Number of rows in an input matrix.
	 * @tparam N
	 *     Number of columns in an input matrix.
	 */
	template < int M, int N >
	struct Svd {
		/**
		 * Tuple of left singular vectors, singular values and right singular
		 * vectors.
		 *
		 * Use `getU`, `getS` and `getV` instead of `std::get`.
		 */
		typedef std::tuple< Matrix< M, M >,
							Matrix< M, N >,
							Matrix< N, N > > USV;

		/** Returns the left singular vectors from a given `USV` tuple. */
		static inline
		const Matrix< M, M >& getU(const std::tuple< Matrix< M, M >,
													 Matrix< M, N >,
													 Matrix< N, N > >& usv)
		{
			return std::get< 0 >(usv);
		}

		/** Returns the singular values from a given `USV` tuple. */
		static inline
		const Matrix< M, N >& getS(const std::tuple< Matrix< M, M >,
													 Matrix< M, N >,
													 Matrix< N, N > >& usv)
		{
			return std::get< 1 >(usv);
		}

		/** Returns the right singular vectors from a given `USV` tuple. */
		static inline
		const Matrix< N, N >& getV(const std::tuple< Matrix< M, M >,
													 Matrix< M, N >,
													 Matrix< N, N > >& usv)
		{
			return std::get< 2 >(usv);
		}

		/**
		 * Decomposes a given matrix into left singular vectors,
		 * singular values and right singular vectors.
		 *
		 * @param m
		 *     `M` x `N` matrix to be decomposed.
		 * @return
		 *     Decomposition of `m`.
		 * @see getU
		 * @see getS
		 * @see getV
		 */
		static USV decomposeUSV(const Matrix< M, N >& m) {
			const int MAX_ITERATIONS = N * 10;
			// allocates matrices
			Matrix< M, M > u = Matrix< M, M >::identity();
			Matrix< M, N > s = m.clone();
			Matrix< N, N > v = Matrix< N, N >::identity();
			// bidiagonalizes a given matrix
			bidiagonalize(u, s, v);
			// repeats Francis iteration
			int n = N;
			for (int i = 0; i < MAX_ITERATIONS; ++i) {
				doFrancis(u, s, v, n);
				// checks the convergence
				if (std::abs(s(n - 2, n - 1) / s(n - 1, n - 1)) < 1.0e-15) {
					--n;
					if (n < 2) {
						std::cout << "iteration: " << i << std::endl;
						break;
					}
				}
			}
			return std::make_tuple(std::move(u), std::move(s), std::move(v));
		}

		/**
		 * Bindiagonalizes a given matrix.
		 *
		 * @param[in,out] u
		 *     Left singular vectors to be upated.
		 * @param[in,out] m
		 *     Matrix to be bidiagonalized.
		 * @param[in,out] v
		 *     Right singular vectors to be updated.
		 */
		static void bidiagonalize(Matrix< M, M >& u,
								  Matrix< M, N >& m,
								  Matrix< N, N >& v)
		{
			for (int i = 0; i < N; ++i) {
				// applies a householder transform to the column vector i
				Reflector< M > rU(m.column(i).slice(i));
				m = rU.applyFromLeftTo(m);
				u = rU.applyFromRightTo(u);
				if (i + 1 < N) {
					// applies a householder transform to the row vector i + 1
					Reflector< N > rV(m.row(i).slice(i + 1));
					m = rV.applyFromRightTo(m);
					v = rV.applyFromRightTo(v);
				}
			}
		}

		/**
		 * Does a single Francis iteration.
		 *
		 * Submatrices other than the top-left `n` x `n` submatrix of `m` are
		 * regarded as already converged.
		 *
		 * @param[in,out] u
		 *     Left singular vectors to be updated.
		 * @param[in,out] m
		 *     Bidiagonalized input matrix to be singular values after
		 *     convergence.
		 * @param[in,out] v
		 *     Right singular vectors to be updated.
		 * @param n
		 *     Size of the submatrix over which the Francis iteration is to be
		 *     performed..
		 */
		static void doFrancis(Matrix< M, M >& u,
							  Matrix< M, N >& m,
							  Matrix< N, N >& v,
							  int n)
		{
			// calculates the shift
			double rho = calculateShift(m, n);
			// applies the first right rotator
			double b1 = m(0, 0);
			double g1 = m(0, 1);
			double mx =
				std::max(std::abs(rho), std::max(std::abs(b1), std::abs(g1)));
			rho /= mx;
			b1 /= mx;
			g1 /= mx;
			Rotator r0(b1 * b1 - rho * rho, b1 * g1);
			m = r0.applyFromRightTo(m, 0);
			v = r0.applyFromRightTo(v, 0);
			// applies the first left rotator
			Rotator r1(m(0, 0), m(1, 0));
			m = r1.applyFromLeftTo(m, 0);
			u = r1.applyFromRightTo(u, 0);
			for (int i = 1; i + 1 < n; ++i) {
				// calculates (i+1)-th right rotator
				Rotator rV(m(i - 1, i), m(i - 1, i + 1));
				m = rV.applyFromRightTo(m, i);
				v = rV.applyFromRightTo(v, i);
				// calculates (i+1)-th left rotator
				Rotator rU(m(i, i), m(i + 1, i));
				m = rU.applyFromLeftTo(m, i);
				u = rU.applyFromRightTo(u, i);
			}
		}

		/**
		 * Calculates the shift for a given bidiagonal matrix.
		 *
		 * Submatrices other than top-left `n` x `n` submatrix of `m` are
		 * regarded as already converged.
		 *
		 * @param m
		 *     Bidiagonal matrix from which a shift is to be calculated.
		 * @param n
		 *     Size of the submatrix to be considered. 
		 * @return
		 *     Shift for the top-left `n` x `n` submatrix of `m`.
		 */
		static double calculateShift(const Matrix< M, N >& m, int n) {
			double b1 = m(n - 2, n - 2);
			double b2 = m(n - 1, n - 1);
			double g1 = m(n - 2, n - 1);
			// solves lambda^4 - d*lambda^2 + e = 0
			// where
			//  d = b1^2 + b2^2 + g1^2
			//  e = b1^2 * b2^2
			// chooses lambda rho closest to b2
			double rho;
			double d = b1 * b1 + b2 * b2 + g1 * g1;
			double e = b1 * b1 * b2 * b2;
			// lambda^2 = (d +- sqrt(d^2 - 4e)) / 2
			// so, f = d^2 - 4e must be positive
			double f = d * d - 4 * e;
			if (f >= 0) {
				f = sqrt(f);
				// lambda = +-sqrt(d +- f)  (d >= 0, f >= 0)
				// if d > f, both d+f and d-f have real square roots
				// otherwise considers only d+f
				if (d > f) {
					// lets l1 > l2
					double l1 = sqrt((d + f) * 0.5);
					double l2 = sqrt((d - f) * 0.5);
					// if b2 >= 0, chooses a positive shift
					// otherwise chooses a negative shift
					if (b2 >= 0) {
						if (std::abs(b2 - l1) < std::abs(b2 - l2)) {
							rho = l1;
						} else {
							rho = l2;
						}
					} else {
						if (std::abs(b2 + l1) < std::abs(b2 + l2)) {
							rho = -l1;
						} else {
							rho = -l2;
						}
					}
				} else {
					double l1 = sqrt((d + f) * 0.5);
					if (std::abs(b2 - l1) <= std::abs(b2 + l1)) {
						rho = l1;
					} else {
						rho = -l1;
					}
				}
			} else {
				// no solution. chooses b2 as the shift
				rho = b2;
			}
			return rho;
		}
	};

}

#endif
