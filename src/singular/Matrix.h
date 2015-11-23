#ifndef _SINGULAR_MATRIX_H
#define _SINGULAR_MATRIX_H

#include "singular/Vector.h"

#include <algorithm>
#include <cstring>
#include <iostream>

namespace singular {

	/**
	 * Fixed sized matrix.
	 *
	 * @tparam M
	 *     Number of rows.
	 * @tparam N
	 *     Number of columns.
	 */
	template < int M, int N >
	class Matrix {
	private:
		/**
		 * Memory block for this matrix.
		 *
		 * Element at the row `i` and column `j` is given by
		 * `pBlock[i * N + j]`.
		 */
		double* pBlock;

		// transposed matrix is a friend
		friend class Matrix< N, M >;
	public:
		/** Initializes a zero matrix. */
		Matrix() {
			this->pBlock = new double[M * N];
			std::fill(this->pBlock, this->pBlock + (M * N), 0.0);
		}

		/** Simple copy is not allowed. */
		Matrix(const Matrix< M, N >& copyee) = delete;

		/**
		 * Steals the memory block of a given matrix.
		 *
		 * @param[in,out] copyee
		 *     Matrix to be copied. No longer valid after this call.
		 */
		Matrix(Matrix< M, N >&& copyee) : pBlock(copyee.pBlock) {
			copyee.pBlock = 0;
		}

		/** Releases the allocated block. */
		~Matrix() {
			this->release();
		}

		/**
		 * Steals the memory block of a given matrix.
		 *
		 * @param[in,out] copyee
		 *     Matrix to be copied. No longer valid after this call.
		 * @return
		 *     Reference to this matrix.
		 */
		Matrix< M, N >& operator =(Matrix< M, N >&& copyee) {
			this->release();
			this->pBlock = copyee.pBlock;
			copyee.pBlock = 0;
			return *this;
		}

		/** Simple copy is not allowed. */
		Matrix< M, N >& operator =(const Matrix< M, N >& copyee) = delete;

		/**
		 * Creates a clone of this matrix.
		 *
		 * A clone has the same contents of this matrix but an independent
		 * memory block from this matrix.
		 *
		 * @return
		 *     Clone of this matrix.
		 */
		Matrix< M, N > clone() const {
			double* pBlock = new double[M * N];
			std::copy(this->pBlock, this->pBlock + M * N, pBlock);
			return Matrix< M, N >(pBlock);
		}

		/**
		 * Creates an identity matrix.
		 *
		 * This function cannot be compiled if `M != N`.
		 *
		 * @return
		 *     Identity matrix.
		 */
		static Matrix< M, N > identity() {
			Matrix< M, M > eye;
			double* pDst = eye.pBlock;
			for (int i = 0; i < M; ++i) {
				*pDst = 1;
				pDst += N + 1;
			}
			return eye;  // compile error if M != N*/
		}

		/**
		 * Creates a matrix filled with given values.
		 *
		 * The value at the ith row and jth column is taken from
		 * `values[M * i + j]`.
		 *
		 * The behavior is undefined if `values` has less than `M * N` elements.
		 *
		 * @param values
		 *     Values to fill the matrix.
		 * @return
		 *     Matrix filled with `values`.
		 */
		static Matrix< M, N > filledWith(const double values[]) {
			double* pBlock = new double[M * N];
			memcpy(pBlock, values, sizeof(double) * M * N);
			return Matrix< M, N >(pBlock);
		}

		/** Returns the value at a given row and column. */
		inline double operator ()(int i, int j) const {
			return this->pBlock[i * N + j];
		}

		/** Returns the value at a given row and column. */
		inline double& operator ()(int i, int j) {
			return this->pBlock[i * N + j];
		}

		/**
		 * Returns a given row in this matrix as a modifiable vector.
		 *
		 * @param i
		 *     Index of the row to be obtained.
		 * @return
		 *     ith row as a vector.
		 *     Changes on this vector are reflected to this matrix.
		 */
		Vector< double > row(int i) {
			return Vector< double >(this->pBlock + i * N, N, 1);
		}

		/**
		 * Returns a given row in this matrix as an unmodifiable vector.
		 *
		 * @param i
		 *     Index of the row to be obtained.
		 * @return
		 *     ith row as a vector.
		 */
		Vector< const double > row(int i) const {
			return Vector< const double >(this->pBlock + i * N, N, 1);
		}

		/**
		 * Returns a given column in this matrix as a modifiable vector.
		 *
		 * @param j
		 *     Index of the column to be obtained.
		 * @return
		 *     jth column as a vector.
		 *     Changes on this vector are reflected to this matrix.
		 */
		Vector< double > column(int j) {
			return Vector< double >(this->pBlock + j, M, N);
		}

		/**
		 * Returns a given column in this matrix as an modifiable vector.
		 *
		 * @param j
		 *     Index of the column to be obtained.
		 * @return
		 *     jth column as a vector.
		 */
		Vector< const double > column(int j) const {
			return Vector< const double >(this->pBlock + j, M, N);
		}

		/**
		 * Fills this matrix with given values.
		 *
		 * The value at the ith row and jth column is taken from
		 * `values[i * N + j]`.
		 *
		 * The behavior is undefined if `values` has less than `M * N` elements.
		 *
		 * @param values
		 *     Values to fill this matrix.
		 */
		Matrix< M, N >& fill(const double values[]) {
			std::copy(values, values + M * N, this->pBlock);
			return *this;
		}

		// Defined outside
		template < int M2, int N2, int L >
		friend Matrix< M2, L > operator *(const Matrix< M2, N2 >& lhs,
										  const Matrix< N2, L >& rhs);

		// Defined outside
		template <
			int M2, int N2, int L, template < int, int > class MatrixLike >
		friend Matrix< M2, L > operator *(const Matrix< M2, N2 >& lhs,
										  const MatrixLike< N2, L >& rhs);

		// Defined outside
		template <
			int M2, int N2, int L, template < int, int > class MatrixLike >
		friend Matrix< M2, L > operator *(const MatrixLike< M2, N2 >& lhs,
										  const Matrix< N2, L >& rhs);

		// Defined outside
		template <
			int M2, int N2, int L,
			template < int, int > class MatrixLike1,
			template < int, int > class MatrixLike2 >
		friend Matrix< M2, L > operator *(const MatrixLike1< M2, N2 >& lhs,
										  const MatrixLike2< N2, L >& rhs);

		/**
		 * Returns the transposition of this matrix.
		 *
		 * @return
		 *     Transposition of this matrix.
		 */
		Matrix< N, M > transpose() const {
			double* pBlock = new double[M * N];
			const double* pSrc = this->pBlock;
			for (int i = 0; i < M; ++i){
				double* pDst = pBlock + i;
				for (int j = 0; j < N; ++j) {
					*pDst = *pSrc;
					++pSrc;
					pDst += M;
				}
			}
			return Matrix< N, M >(pBlock);
		}

		/**
		 * Shuffles rows in this matrix.
		 *
		 * Equivalent to multiplying the following permutation matrix \f$P\f$
		 * from the left of this matrix.
		 *
		 * \f[
		 * P_{ij} =
		 *   \left{
		 *     \begin{array}{ll}
		 *       1 & (order[i] = j) \\
		 *       0 & (order[i] \neq j)
		 *     \end{array}
		 *   \right.
		 * \f]
		 *
		 * @param order
		 *     New order of rows.
		 *     Must have at least M elements.
		 * @return
		 *     Matrix shuffled in the given order.
		 */
		Matrix< M, N > shuffleRows(const int order[]) const {
			double* pBlock = new double[M * N];
			double* pDst = pBlock;
			for (int i = 0; i < M; ++i) {
				double* pSrc = this->pBlock + order[i] * N;
				std::copy(pSrc, pSrc + N, pDst);
				pDst += N;
			}
			return Matrix< M, N >(pBlock);
		}

		/**
		 * Shuffles columns in this matrix.
		 *
		 * Equivalent to multiplying the following permutation matrix \f$P\f$
		 * from the right of this matrix.
		 *
		 * \f[
		 * P_{ij} =
		 *   \left{
		 *     \begin{array}{ll}
		 *       1 & (i = order[j]) \\
		 *       0 & (i \neq order[j])
		 *     \end{array}
		 *   \right.
		 * \f]
		 *
		 * @param order
		 *     New order of columns.
		 *     Must have at least N elements.
		 * @return
		 *     Matrix shuffled in the given order.
		 */
		Matrix< M, N > shuffleColumns(const int order[]) const {
			double* pBlock = new double[M * N];
			for (int j = 0; j < N; ++j) {
				double* pDst = pBlock + j;
				double* pSrc = this->pBlock + order[j];
				for (int i = 0; i < M; ++i) {
					*pDst = *pSrc;
					pSrc += N;
					pDst += N;
				}
			}
			return Matrix< M, N >(pBlock);
		}
	private:
		/**
		 * Initializes with a given memory block.
		 *
		 * @param pBlock
		 *     Memory block of the new matrix.
		 *     Must have at least M * N elements.
		 */
		Matrix(double* pBlock) : pBlock(pBlock) {}

		/**
		 * Releases the memory block of this matrix.
		 *
		 * Has no effect if the memory block has already been released.
		 */
		inline void release() {
			delete[] this->pBlock;
			this->pBlock = 0;
		}
	};

	/**
	 * Multiplies given two matrices.
	 *
	 * @tparam M
	 *     Number of rows in the left-hand-side matrix.
	 * @tparam N
	 *     Number of columns in the left-hand-side matrix.
	 *     Number of rows in the right-hand-side matrix as well.
	 * @tparam L
	 *     Number of columns in the right-hand-side matrix.
	 * @param lhs
	 *     Left-hand side of the multiplication.
	 * @param rhs
	 *     Right-hand side of the multiplication.
	 * @return
	 *     Product of `lhs` and `rhs`.
	 */
	template < int M, int N, int L >
	Matrix< M, L > operator *(const Matrix< M, N >& lhs,
							  const Matrix< N, L >& rhs)
	{
		double* pBlock = new double[M * L];
		double* pDst = pBlock;
		for (int i = 0; i < M; ++i) {
			for (int l = 0; l < L; ++l) {
				double* pL = lhs.pBlock + i * N;
				double* pR = rhs.pBlock + l;
				double x = 0.0;
				for (int j = 0; j < N; ++j) {
					x += *pL * *pR;
					++pL;
					pR += L;
				}
				*pDst = x;
				++pDst;
			}
		}
		return Matrix< M, L >(pBlock);
	}

	/**
	 * Multiplies given two matrices.
	 *
	 * `MatrixLike` must overload the function-call operator `()` which
	 * takes a row index `i` and column index `j` and returns the element at
	 * the ith row and jth column.
	 * The function prototype should look like the following,
	 *  - `double operator ()(int i, int j) const`
	 *
	 * @tparam M
	 *     Number of rows in the left-hand-side matrix.
	 * @tparam N
	 *     Number of columns in the left-hand-side matrix.
	 *     Number of rows in the right-hand-side matrix.
	 * @tparam L
	 *     Number of columns in the right-hand-side matrix.
	 * @tparam MatrixLike
	 *     Type of the right-hand-side matrix.
	 * @param lhs
	 *     Left-hand-side of the multiplication.
	 * @param rhs
	 *     Right-hand-side of the multiplication.
	 * @return
	 *     Product of `lhs` and `rhs`.
	 */
	template < int M, int N, int L, template < int, int > class MatrixLike >
	Matrix< M, L > operator *(const Matrix< M, N >& lhs,
							  const MatrixLike< N, L >& rhs)
	{
		double* pBlock = new double[M * L];
		double* pDst = pBlock;
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < L; ++j) {
				const double* pL = lhs.pBlock + i * N;
				double x = 0.0;
				for (int k = 0; k < N; ++k) {
					x += *pL * rhs(k, j);
					++pL;
				}
				*pDst = x;
				++pDst;
			}
		}
		return Matrix< M, L >(pBlock);
	}

	/**
	 * Multiplies given two matrices.
	 *
	 * `MatrixLike` must overload the function-call operator `()` which
	 * takes a row index `i` and column index `j` and returns the element at
	 * the ith row and jth column.
	 * The function prototype should look like the following,
	 *  - `double operator ()(int i, int j) const`
	 *
	 * @tparam M
	 *     Number of rows in the left-hand-side matrix.
	 * @tparam N
	 *     Number of columns in the left-hand-side matrix.
	 *     Number of rows in the right-hand-side matrix.
	 * @tparam L
	 *     Number of columns in the right-hand-side matrix.
	 * @tparam MatrixLike
	 *     Type of the left-hand-side matrix.
	 * @param lhs
	 *     Left-hand-side of the multiplication.
	 * @param rhs
	 *     Right-hand-side of the multiplication.
	 * @return
	 *     Product of `lhs` and `rhs`.
	 */
	template < int M, int N, int L, template < int, int > class MatrixLike >
	Matrix< M, L > operator *(const MatrixLike< M, N >& lhs,
							  const Matrix< N, L >& rhs)
	{
		double* pBlock = new double[M * L];
		double* pDst = pBlock;
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < L; ++j) {
				const double* pR = rhs.pBlock + j;
				double x = 0.0;
				for (int k = 0; k < N; ++k) {
					x += lhs(i, k) * *pR;
					pR += L;
				}
				*pDst = x;
				++pDst;
			}
		}
		return Matrix< M, L >(pBlock);
	}

	/**
	 * Multiplies given two matrices.
	 *
	 * Both of `MatrixLike1` and `MatrixLike2` must overload the
	 * function-call operator `()` which takes a row index `i` and column
	 * index `j` and returns the element at ith row and jth column.
	 * The function prototype should look like the following,
	 *  - `double operator ()(int i, int j) const`
	 *
	 * @tparam M
	 *     Number of rows in the left-hand-side matrix.
	 * @tparam N
	 *     Number of columns in the left-hand-side matrix.
	 *     Number of rows in the right-hand-side matrix.
	 * @tparam L
	 *     Number of columns in the right-hand-side matrix.
	 * @tparam MatrixLike1
	 *     Type of the left-hand-side matrix.
	 * @tparam MatrixLike2
	 *     Type of the right-hand-side matrix.
	 * @param lhs
	 *     Left-hand-side of the multiplication.
	 * @param rhs
	 *     Right-hand-side of the multiplication.
	 * @return
	 *     Product of `lhs` and `rhs`.
	 */
	template <
		int M, int N, int L,
		template < int, int > class MatrixLike1,
		template < int, int > class MatrixLike2 >
	Matrix< M, L > operator *(const MatrixLike1< M, N >& lhs,
							  const MatrixLike2< N, L >& rhs)
	{
		double* pBlock = new double[M * L];
		double* pDst = pBlock;
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < L; ++j) {
				double x = 0.0;
				for (int k = 0; k < N; ++k) {
					x += lhs(i, k) * rhs(k, j);
				}
				*pDst = x;
				++pDst;
			}
		}
		return Matrix< M, L >(pBlock);
	}


	/**
	 * Writes a given matrix to a given stream.
	 *
	 * `MatrixLike` should overload the function-call operator `()` which takes
	 * a row index `i` and column index `j` and returns the element at the ith
	 * row and jth column.
	 * The function prototype should look like the following,
	 *  - `double operator ()(int i, int j) const`
	 *
	 * @tparam M
	 *     Number of rows in the matrix to be written.
	 * @tparam N
	 *     Number of columns in the matrix to be written.
	 * @tparam MatrixLike
	 *     Type of the matrix to be written.
	 * @param out
	 *     Output stream where the matrix is to be written.
	 * @param m
	 *     Matrix to be written.
	 * @return
	 *     `out`.
	 */
	template < int M, int N, template < int, int > class MatrixLike >
	std::ostream& operator <<(std::ostream& out, const MatrixLike< M, N >& m) {
		out << '[' << std::endl;
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < N; ++j) {
				out << m(i, j);
				if (j + 1 < N) {
					out << ' ';
				}
			}
			out << std::endl;
		}
		out << ']';
		return out;
	}

}

#endif
