#ifndef _SINGULAR_SINGULAR_REFLECTOR_H
#define _SINGULAR_SINGULAR_REFLECTOR_H

#include "singular/Matrix.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <vector>

namespace singular {

	/**
	 * Reflector.
	 *
	 * A reflector transforms a vector \f$(\mathbf{x} = [x_1 x_2 ... x_N])
	 * into a vector \f$(\mathbf{y} = [-\tau 0 ... 0]).
	 *
	 * A reflector is equivalent to the following matrix.
	 * \f[
	 * \mathbf{H} = \mathbf{I} - \gamma \mathbf{u} \mathbf{u}^T
	 * \f]
	 *
	 * @tparam L
	 *     Size of the transform matrix.
	 */
	template < int L >
	class Reflector {
	private:
		/** U vector. */
		std::vector< double > u;

		/** Gamma. */
		double gamma;
	public:
		/**
		 * Constructs a reflector from a given vector.
		 *
		 * Forms an `L x L` matrix like the following,
		 *
		 * \f[
		 * \left|
		 * \begin{array}{cc}
		 * \mathbf{I} & \mathbf{0} \\
		 * \mathbf{0} & \mathbf{R} \\
		 * \end{array}
		 * \right|
		 * \f]
		 *
		 * ```
		 * | I | 0 |
		 * |---+---|
		 * | 0 | R |
		 * ```
		 *
		 * \f$\mathbf{R}\f$ is an `N x N` reflector created from `v` where
		 * `N = v.size()`.
		 * \f$\mathbf{I}\f$ is an `(L-N)x(L-N)` identity matrix.
		 *
		 * @param v
		 *     Vector from which the reflector is formed.
		 * @throws Exception
		 *     If `v.size() == 0`, or if `v.size() > L`.
		 */
		Reflector(const Vector< const double >& v) {
			const size_t N = v.size();
			// copies the vector
			this->u.reserve(N);
			std::copy(v.begin(), v.end(), std::back_inserter(this->u));
			// normalizes elements by the maximum amplitude
			// to avoid harmful underflow and overflow
			double mx = 0.0;
			for (size_t i = 0; i < N; ++i) {
				mx = std::max(std::abs(this->u[i]), mx);
			}
			if (mx > 0.0) {
				// calculates the normalized norm
				double tau = 0.0;
				for (size_t i = 0; i < N; ++i) {
					double x = this->u[i] / mx;
					this->u[i] = x;
					tau += x * x;
				}
				tau = sqrt(tau);
				// tau's sign should be the same as the first element in `u`
				if (this->u[0] < 0.0) {
					tau = -tau;
				}
				double u0 = this->u[0] + tau;
				this->u[0] = u0;
				std::transform(this->u.begin(), this->u.end(), this->u.begin(),
					[u0] (double& e) { return e / u0; });
				this->gamma = u0 / tau;
			} else {
				// v is a zero vector
				this->gamma = 0.0;
				std::fill(this->u.begin(), this->u.end(), 0.0);
			}
		}

		/**
		 * Applies this reflector to a given matrix from left.
		 *
		 * @tparam N
		 *     Number of columns in the given matrix.
		 * @param m
		 *     Matrix to be transformed.
		 * @return
		 *     Transformed matrix.
		 */
		template < int N >
		Matrix< L, N > applyFromLeftTo(const Matrix< L, N >& m) const {
			// H * m = m - gamma * u * u^T * m
			Matrix< L, N > m2 = m.clone();
			int offset = L - u.size();
			for (int i = 0; i < N; ++i) {
				// caches gamma * u^T * m
				double gUM = 0.0;
				for (int j = 0; j < u.size(); ++j) {
					gUM += this->u[j] * m(j + offset, i);
				}
				gUM *= this->gamma;
				// H * m = m - u * gUM
				for (int j = 0; j < u.size(); ++j) {
					m2(j + offset, i) = m(j + offset, i) - (this->u[j] * gUM);
				}
			}
			return m2;
		}

		/**
		 * Applies this reflector to a given matrix from right.
		 *
		 * @tparam M
		 *     Number of rows in the given matrix.
		 * @param m
		 *     Matrix to be transformed.
		 * @return
		 *     Transformed matrix.
		 */
		template < int M >
		Matrix< M, L > applyFromRightTo(const Matrix< M, L >& m) const {
			// m * H = m - m * gamma * u * u^T
			Matrix< M, L > m2 = m.clone();
			int offset = L - u.size();
			for (int i = 0; i < M; ++i) {
				// caches gamma * m * u
				double gMU = 0.0;
				for (int j = 0; j < u.size(); ++j) {
					gMU += m(i, j + offset) * this->u[j];
				}
				gMU *= this->gamma;
				// m * H = m - gMU * u^T
				for (int j = 0; j < u.size(); ++j) {
					m2(i, j + offset) = m(i, j + offset) - (gMU * this->u[j]);
				}
			}
			return m2;
		}
	private:
		/**
		 * Prints a given reflector to a given output stream.
		 *
		 * @param out
		 *     Output stream to which `reflector` is to be printed.
		 * @param reflector
		 *     Reflector to be printed.
		 * @return
		 *     `out`.
		 */
		friend std::ostream& operator <<(std::ostream& out,
										 const Reflector< L >& reflector)
		{
			out << "gamma: " << reflector.gamma;
			out << ", u: [";
			for (size_t i = 0; i < reflector.u.size(); ++i) {
				out << reflector.u[i];
				if (i + 1 < reflector.u.size()) {
					out << " ";
				}
			}
			out << "]";
			return out;
		}
	};

}

#endif
