#ifndef _SINGULAR_SVD_H
#define _SINGULAR_SVD_H

#include "singular/Matrix.h"
#include "singular/Reflector.h"
#include "singular/Rotator.h"

#include <tuple>

namespace singular {

	/**
	 * Performs singular value decomposition of a given matrix.
	 *
	 * `M`: Number of rows.
	 *
	 * `N`: Number of columns.
	 *
	 * @tparam M
	 *     Number of rows of the input matrix.
	 * @tparam N
	 *     Number of columns of the input matrix.
	 * @param m
	 *     `M x N` matrix to be decomposed.
	 * @return
	 *     Decomposition of `m`.
	 */
	template < int M, int N >
	std::tuple< Matrix< M, M >, Matrix< M, N >, Matrix< N, N > >
		svdUSV(const Matrix< M, N >& m)
	{
		// bidiagonalizes a given matrix
		std::tuple< Matrix< M, M >,
					Matrix< M, N >,
					Matrix< N, N > > usv = bidiagonalize(m);
		// does a single Francis iteration
		usv = doFrancis(usv);
		return usv;
		/*
		return std::make_tuple(std::move(u),
							   Matrix< M, N >(),
							   std::move(v));*/
	}

	/**
	 * Bindiagonalizes a given matrix.
	 *
	 * @tparam M
	 *     Number of rows of the input matrix.
	 * @tparam N
	 *     Number of columns of the input matrix.
	 * @param m
	 *     `M x N` matrix to be bidiagonalized.
	 * @return
	 *     Bidiagonalized matrix.
	 */
	template < int M, int N >
	std::tuple< Matrix< M, M >, Matrix< M, N >, Matrix< N, N > >
		bidiagonalize(const Matrix< M, N >& m)
	{
		Matrix< M, N > s = m.clone();
		Matrix< M, M > u = Matrix< M, M >::identity();
		Matrix< N, N > v = Matrix< N, N >::identity();
		for (int i = 0; i < N; ++i) {
			// applies a householder transform to the column vector i
			Reflector< M > colReflector(s.column(i).slice(i));
			s = colReflector.applyFromLeftTo(s);
			u = colReflector.applyFromLeftTo(u);
			if (i + 1 < N) {
				// applies a householder transform to the row vector i + 1
				Reflector< N > rowReflector(s.row(i).slice(i + 1));
				s = rowReflector.applyFromRightTo(s);
				v = rowReflector.applyFromRightTo(v);
			}
		}
		return std::make_tuple(std::move(u), std::move(s), std::move(v));
	}

	/**
	 * Does a single Francis iteration.
	 *
	 * @tparam M
	 *     Number of rows in the input matrix.
	 * @tparam N
	 *     Number of columns in the input matrix.
	 * @param ubv
	 *     Tuple of U, B, V matrices.
	 *     B must be a bidiagonalized input matrix.
	 * @return
	 *     Tuple of updated U, B(~S), V matrices.
	 */
	template < int M, int N >
	std::tuple< Matrix< M, M >, Matrix< M, N >, Matrix< N, N > >
		doFrancis(std::tuple< Matrix< M, M >,
							  Matrix< M, N >,
							  Matrix< N, N > >& ubv)
	{
		Matrix< M, M > u = std::get< 0 >(ubv).clone();
		Matrix< M, N > b = std::get< 1 >(ubv).clone();
		Matrix< N, N > v = std::get< 2 >(ubv).clone();
		// calculates the shift
		double lo = getShift(b);
		// applies the first right rotator
		double b1 = b(0, 0);
		double g1 = b(0, 1);
		double mx =
			std::max(std::abs(lo), std::max(std::abs(b1), std::abs(g1)));
		lo /= mx;
		b1 /= mx;
		g1 /= mx;
		Rotator r0(b1 * b1 - lo * lo, b1 * g1);
		b = r0.applyFromRightTo(b, 0);
		v = r0.applyFromRightTo(v, 0);
		// applies the first left rotator
		Rotator r1(b(0, 0), b(1, 0));
		b = r1.applyFromLeftTo(b, 0);
		u = r1.applyFromLeftTo(u, 0);
		for (int i = 1; i + 1 < N; ++i) {
			// calculates (i+1)-th right rotator
			Rotator rV(b(i - 1, i), b(i - 1, i + 1));
			b = rV.applyFromRightTo(b, i);
			v = rV.applyFromRightTo(v, i);
			// calculates (i+1)-th left rotator
			Rotator rU(b(i, i), b(i + 1, i));
			b = rU.applyFromLeftTo(b, i);
			u = rU.applyFromLeftTo(u, i);
		}
		return std::make_tuple(std::move(u), std::move(b), std::move(v));
	}

	/**
	 * Returns the shift for a given bidiagonal matrix.
	 *
	 * @tparam M
	 *     Number of rows in the given bidiagonal matrix.
	 * @tparam N
	 *     Number of columns in the given bidiagonal matrix.
	 * @param b
	 *     Bidiagonal matrix from which a shift is obtained.
	 * @return
	 *     Shift for `b`.
	 */
	template < int M, int N >
	double getShift(const Matrix< M, N >& b) {
		return b(N - 1, N - 1);
	}

}

#endif
