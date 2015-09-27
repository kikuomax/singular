#ifndef _SINGULAR_ROTATOR_H
#define _SINGULAR_ROTATOR_H

#include "singular/Matrix.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace singular {

	/**
	 * Rotator.
	 */
	class Rotator {
	private:
		/** Cosine value of this rotator. */
		double cs;

		/** Sine value of this rotator. */
		double sn;
	public:
		/**
		 * Builds a rotator from a given two-element vector.
		 *
		 * Builds a `2x2` rotator `Q` such that,
		 *
		 * ```
		 * Q^T *| x1 | = | * |
		 *      | x2 |   | 0 |
		 * ```
		 *
		 * @param x1
		 *     First element in the vector.
		 * @param x2
		 *     Second element in the vector.
		 * @param offset
		 *     Row index of `x1`.
		 */
		Rotator(double x1, double x2) {
			// normalizes by the maximum magnitude
			// to avoid harmful underflow and overflow
			double mx = std::max(std::abs(x1), std::abs(x2));
			x1 /= mx;
			x2 /= mx;
			double norm = sqrt(x1 * x1 + x2 * x2);
			this->cs = x1 / norm;
			this->sn = x2 / norm;
		}

		/**
		 * Applies this rotator from the left hand side of a given matrix.
		 *
		 * This rotator can be viewed as the following MxM matrix,
		 *
		 * ```
		 * | I1        |
		 * |    Q^T    |
		 * |        I2 |
		 *
		 * I1: k x k identity matrix
		 * I2: (M-(k+2))x(M-(k+2)) identity matrix
		 * ```
		 *
		 * Undefined if `M < k + 2`.
		 * 
		 * @tparam M
		 *     Number of the rows in the given matrix.
		 * @tparam N
		 *     Number of the columns in the given matrix.
		 * @param rhs
		 *     Matrix to be rotated.
		 * @param k
		 *     Row and column index where this rotator is applied.
		 * @return
		 *     Result of this rotation.
		 */
		template < int M, int N >
		Matrix< M, N > applyFromLeftTo(const Matrix< M, N >& rhs, int k) {
			Matrix< M, N > m = rhs.clone();
			for (int i = 0; i < N; ++i) {
				double x1 = rhs(k, i);
				double x2 = rhs(k + 1, i);
				m(k, i) = this->cs * x1 + this->sn * x2;
				m(k + 1, i) = -this->sn * x1 + this->cs * x2;
			}
			return m;
		}

		/**
		 * Applies this rotator from the right hand side of a given matrix.
		 *
		 * This rotator can be viewed as the following NxN matrix,
		 *
		 * ```
		 * | I1      |
		 * |    Q    |
		 * |      I2 |
		 *
		 * I1: k x k identity matrix
		 * I2: (N-(k+2))x(N-(k+2)) identity matrix
		 * ```
		 *
		 * Undefined if `N < k + 2`.
		 *
		 * @tparam M
		 *     Number of the rows in the given matrix.
		 * @tparam N
		 *     Number of the columns in the given matrix.
		 * @param lhs
		 *     Matrix to be rotated.
		 * @param k
		 *     Row and column index where this rotator is applied.
		 * @return
		 *     Result of this rotation.
		 */
		template < int M, int N >
		Matrix< M, N > applyFromRightTo(const Matrix< M, N >& lhs, int k) {
			Matrix< M, N > m = lhs.clone();
			for (int i = 0; i < M; ++i) {
				double x1 = lhs(i, k);
				double x2 = lhs(i, k + 1);
				m(i, k) = x1 * this->cs + x2 * this->sn;
				m(i, k + 1) = x1 * (-this->sn) + x2 * this->cs;
			}
			return m;
		}
	private:
	};

}

#endif
