#ifndef _SINGULAR_MATRIX_H
#define _SINGULAR_MATRIX_H

#include "singular/Vector.h"

#include <algorithm>
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
		/** Memory block for this matrix. */
		double* pBlock;
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
		virtual ~Matrix() {
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
			Matrix cloneM;
			memcpy(cloneM.pBlock, this->pBlock, sizeof(double) * M * N);
			return cloneM;
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
			for (int i = 0; i < M; ++i) {
				eye(i, i) = 1;
			}
			return eye;  // compile error if M != N
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
		 * Returns a given row in this matrix as a vector.
		 *
		 * @param i
		 *     Index of the row to be obtained.
		 * @return
		 *     i-th row as a vector.
		 */
		Vector< const double > row(int i) const {
			return Vector< const double >(this->pBlock + i * N, N, 1);
		}

		/**
		 * Returns a given column in this matrix as a vector.
		 *
		 * @param i
		 *     Index of the column to be obtained.
		 * @return
		 *     i-th column as a vector.
		 */
		Vector< const double > column(int i) const {
			return Vector< const double >(this->pBlock + i, M, N);
		}

		/**
		 * Fills this matrix with given values.
		 *
		 * The value at row `i` and column `j` will be taken from
		 * `values[i * N + j]`.
		 *
		 * @param values
		 *     Array of values to fill this matrix.
		 *     Must have at least `M * N` valid elements.
		 */
		Matrix< M, N >& fill(const double values[]) {
			for (int i = 0; i < M; ++i) {
				for (int j = 0; j < N; ++j) {
					this->pBlock[i * N + j] = values[i * N + j];
				}
			}
			return *this;
		}
	private:
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

	/** Writes a given matrix to a given stream. */
	template < int M, int N >
	std::ostream& operator <<(std::ostream& out, const Matrix< M, N >& m) {
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
