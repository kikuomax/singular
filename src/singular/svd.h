#ifndef _SINGULAR_SINGULAR_SVD_H
#define _SINGULAR_SINGULAR_SVD_H

#include "singular/Matrix.h"
#include "singular/Reflector.h"

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
		/*
		// allocates U, V matrices as identity matrices
		Matrix< M, M > u = Matrix< M, M >::identity();
		Matrix< N, N > v = Matrix< N, N >::identity();*/
		// bidiagonalizes a given matrix
		std::tuple< Matrix< M, M >,
					Matrix< M, N >,
					Matrix< N, N > > usv = bidiagonalize(m);
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
		Matrix< M, N > m2 = m.clone();
		Matrix< M, M > u = Matrix< M, M >::identity();
		Matrix< N, N > v = Matrix< N, N >::identity();
		for (int i = 0; i + 1 < N; ++i) {
			// applies a householder transform to the column vector i
			Reflector< M > colReflector(m.column(i).slice(i));
			m2 = colReflector.applyFromLeftTo(m2);
			u = colReflector.applyFromLeftTo(u);
			if (i + 2 <= N) {
				// applies a householder transform to the row vector i + 1
				Reflector< N > rowReflector(m.row(i + 1).slice(i + 1));
				m2 = rowReflector.applyFromRightTo(m2);
				v = rowReflector.applyFromRightTo(v);
			}
		}
		return std::make_tuple(std::move(u), std::move(m2), std::move(v));
	}

}

#endif
