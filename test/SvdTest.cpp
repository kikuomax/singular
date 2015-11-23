#include "singular/Svd.h"

#include "gtest/gtest.h"

/** Fixture for SVD on a 5x4 matrix. */
class SvdOn5x4MatrixTest : public ::testing::Test {
protected:
	/** Number of rows in the input matrix. */
	static const int M = 5;

	/** Number of columns in the input matrix. */
	static const int N = 4;

	/** Input matrix. */
	singular::Matrix< M, N > m;

	/** Results of SVD. */
	singular::Svd< M, N >::USV usv;

	/** Builds a matrix and performs SVD on it. */
	virtual void SetUp() {
		const double DATA[] = {
			1.0, 2.0, 3.0, 4.0,
			5.0, 6.0, 7.0, 8.0,
			4.0, 8.0, 3.0, 5.0,
			6.0, 7.0, 2.0, 1.0,
			9.0, 1.0, 3.0, 6.0
		};
		this->m.fill(DATA);
		this->usv = singular::Svd< M, N >::decomposeUSV(m);
	}
};

TEST_F(SvdOn5x4MatrixTest, Left_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	// U*U^T = I
	singular::Matrix< M, M > eye =
		singular::Svd< M, N >::getU(this->usv)
		* singular::Svd< M, N >::getU(this->usv).transpose();
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < M; ++j) {
			if (i == j) {
				EXPECT_NEAR(1.0, eye(i, j), ROUNDED_ERROR);
			} else {
				EXPECT_NEAR(0.0, eye(i, j), ROUNDED_ERROR);
			}
		}
	}
}

TEST_F(SvdOn5x4MatrixTest, Right_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	// V*V^T = I
	singular::Matrix< N, N > eye =
		singular::Svd< M, N >::getV(this->usv)
		* singular::Svd< M, N >::getV(this->usv).transpose();
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (i == j) {
				EXPECT_NEAR(1.0, eye(i, j), ROUNDED_ERROR);
			} else {
				EXPECT_NEAR(0.0, eye(i, j), ROUNDED_ERROR);
			}
		}
	}
}

TEST_F(SvdOn5x4MatrixTest, Multiplication_of_USV_should_be_input_matrix) {
	const double ROUNDED_ERROR = 1.0e-13;
	// m = U*S*V^T
	singular::Matrix< M, N > m2 =
		singular::Svd< M, N >::getU(this->usv)
		* singular::Svd< M, N >::getS(this->usv)
		* singular::Svd< M, N >::getV(this->usv).transpose();
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			EXPECT_NEAR(this->m(i, j), m2(i, j), ROUNDED_ERROR);
		}
	}
}

TEST_F(SvdOn5x4MatrixTest, 4_positive_singular_values_should_be_produced) {
	const double ROUNDED_ERROR = 1.0e-14;
	const singular::DiagonalMatrix< M, N >& s =
		singular::Svd< M, N >::getS(this->usv);
	EXPECT_NEAR(21.3113428837071, s(0, 0), ROUNDED_ERROR * 10);
	EXPECT_NEAR(6.71730295404777, s(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(5.77467474261999, s(2, 2), ROUNDED_ERROR);
	EXPECT_NEAR(1.53545990945876, s(3, 3), ROUNDED_ERROR);
}

/** Fixture for SVD on a 4x5 matrix. */
class SvdOn4x5MatrixTest : public ::testing::Test {
protected:
	/** Number of rows in the input matrix. */
	static const int M = 4;

	/** Number of columns in the input matrix. */
	static const int N = 5;

	/** Input matrix. */
	singular::Matrix< M, N > m;

	/** Results of SVD. */
	singular::Svd< M, N >::USV usv;

	/** Builds an input matrix and performs SVD on it. */
	virtual void SetUp() {
		const double DATA[] = {
			3.5, -0.4, 2.7, 1.5, 5.0,
			-2.0, 9.2, 1.1, 0.5, 3.8,
			4.9, 5.5, 4.7, -2.9, 6.0,
			8.2, 1.3, 5.4, 2.6, -1.0
		};
		this->m.fill(DATA);
		this->usv = singular::Svd< M, N >::decomposeUSV(this->m);
	}
};

TEST_F(SvdOn4x5MatrixTest, Left_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	// U*U^T = I
	const singular::Matrix< M, M > eye =
		singular::Svd< M, N >::getU(this->usv)
		* singular::Svd< M, N >::getU(this->usv).transpose();
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < M; ++j) {
			if (i == j) {
				EXPECT_NEAR(1.0, eye(i, j), ROUNDED_ERROR);
			} else {
				EXPECT_NEAR(0.0, eye(i, j), ROUNDED_ERROR);
			}
		}
	}
}

TEST_F(SvdOn4x5MatrixTest, Right_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	// V*V^T = I
	const singular::Matrix< N, N > eye =
		singular::Svd< M, N >::getV(this->usv)
		* singular::Svd< M, N >::getV(this->usv).transpose();
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (i == j) {
				EXPECT_NEAR(1.0, eye(i, j), ROUNDED_ERROR);
			} else {
				EXPECT_NEAR(0.0, eye(i, j), ROUNDED_ERROR);
			}
		}
	}
}

TEST_F(SvdOn4x5MatrixTest, Multiplication_of_USV_should_be_input_matrix) {
	const double ROUNDED_ERROR = 1.0e-13;
	// m = U*S*V^T
	const singular::Matrix< M, N > m2 =
		singular::Svd< M, N >::getU(this->usv)
		* singular::Svd< M, N >::getS(this->usv)
		* singular::Svd< M, N >::getV(this->usv).transpose();
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			EXPECT_NEAR(this->m(i, j), m2(i, j), ROUNDED_ERROR);
		}
	}
}

TEST_F(SvdOn4x5MatrixTest, 4_positive_singular_values_should_be_produced) {
	// reference values were calculated at the following site
	// http://www.bluebit.gr/matrix-calculator/
	const double ROUNDED_ERROR = 1.0e-14;
	const singular::DiagonalMatrix< M, N >& s =
		singular::Svd< M, N >::getS(this->usv);
	EXPECT_NEAR(15.0341927179405, s(0, 0), ROUNDED_ERROR * 10);
	EXPECT_NEAR(10.5443636529196, s(1, 1), ROUNDED_ERROR * 10);
	EXPECT_NEAR(5.37588907505156, s(2, 2), ROUNDED_ERROR);
	EXPECT_NEAR(3.46255124547698, s(3, 3), ROUNDED_ERROR);
}

/** Fixture for SVD on 3x3 but rank 2 matrix. */
class SvdOn3x3Rank2MatrixTest : public ::testing::Test {
protected:
	/** Number of rows in the input matrix. */
	static const int M = 3;

	/** Number of columns in the input matrix. */
	static const int N = 3;

	/** Input matrix. */
	singular::Matrix< M, N > m;

	/** Results of SVD. */
	singular::Svd< M, N >::USV usv;

	/** Builds an input matrix and performs SVD. */
	virtual void SetUp() {
		const double DATA[] = {
			1.0,  1.0, 3.0,
			2.0, -5.0, 4.0,
			1.0,  1.0, 3.0
		};
		this->m.fill(DATA);
		this->usv = singular::Svd< M, N >::decomposeUSV(this->m);
	}
};

TEST_F(SvdOn3x3Rank2MatrixTest, Left_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	// U*U^T = I
	singular::Matrix< M, M > eye =
		singular::Svd< M, N >::getU(this->usv)
		* singular::Svd< M, N >::getU(this->usv).transpose();
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < M; ++j) {
			if (i == j) {
				EXPECT_NEAR(1.0, eye(i, j), ROUNDED_ERROR);
			} else {
				EXPECT_NEAR(0.0, eye(i, j), ROUNDED_ERROR);
			}
		}
	}
}

TEST_F(SvdOn3x3Rank2MatrixTest, Right_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	// V*V^T = I
	singular::Matrix< N, N > eye =
		singular::Svd< M, N >::getV(this->usv)
		* singular::Svd< M, N >::getV(this->usv).transpose();
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (i == j) {
				EXPECT_NEAR(1.0, eye(i, j), ROUNDED_ERROR);
			} else {
				EXPECT_NEAR(0.0, eye(i, j), ROUNDED_ERROR);
			}
		}
	}
}

TEST_F(SvdOn3x3Rank2MatrixTest, Multiplication_of_USV_should_be_input_matrix) {
	const double ROUNDED_ERROR = 1.0e-13;
	// A = U*S*V^T
	singular::Matrix< M, N > m2 =
		singular::Svd< M, N >::getU(this->usv)
		* singular::Svd< M, N >::getS(this->usv)
		* singular::Svd< M, N >::getV(this->usv).transpose();
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			EXPECT_NEAR(this->m(i, j), m2(i, j), ROUNDED_ERROR);
		}
	}
}

TEST_F(SvdOn3x3Rank2MatrixTest, 2_positive_and_1_zero_singular_values_should_be_produced) {
	const double ROUNDED_ERROR = 1.0e-14;
	const singular::DiagonalMatrix< M, N >& s =
		singular::Svd< M, N >::getS(this->usv);
	EXPECT_NEAR(7.11714246017378, s(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(4.04305369758943, s(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, s(2, 2), ROUNDED_ERROR);
}

/** Fixture for SVD on a 7x5 sparse matrix. */
class SvdOn7x5SparseMatrixTest : public ::testing::Test {
protected:
	/** Number of rows in the input matrix. */
	static const int M = 7;

	/** Number of columns in the input matrix. */
	static const int N = 5;

	/** Input matrix. */
	singular::Matrix< M, N > m;

	/** Results of SVD. */
	singular::Svd< M, N >::USV usv;

	/** Builds an input matrix and performs SVD on it. */
	virtual void SetUp() {
		const double DATA[] = {
			1.5, 0.0, 2.0, 1.0, 0.0,
			0.0, -4.2, 0.0, 0.0, 7.6,
			0.0, 8.7, 0.0, 0.0, -2.4,
			7.8, 0.0, 3.7, 4.5, 0.0,
			4.1, 0.0, -1.3, -6.1, 0.0,
			3.2, 0.0, 5.1, 5.5, 0.0,
			2.1, 0.0, -8.9, 9.1, 0.0
		};
		this->m.fill(DATA);
		this->usv = singular::Svd< M, N >::decomposeUSV(this->m);
	}
};

TEST_F(SvdOn7x5SparseMatrixTest, Left_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	singular::Matrix< M, M > eye =
		singular::Svd< M, N >::getU(this->usv)
		* singular::Svd< M, N >::getU(this->usv).transpose();
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < M; ++j) {
			if (i == j) {
				EXPECT_NEAR(1.0, eye(i, j), ROUNDED_ERROR);
			} else {
				EXPECT_NEAR(0.0, eye(i, j), ROUNDED_ERROR);
			}
		}
	}
}

TEST_F(SvdOn7x5SparseMatrixTest, Right_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	singular::Matrix< N, N > eye =
		singular::Svd< M, N >::getV(this->usv)
		* singular::Svd< M, N >::getV(this->usv).transpose();
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (i == j) {
				EXPECT_NEAR(1.0, eye(i, j), ROUNDED_ERROR);
			} else {
				EXPECT_NEAR(0.0, eye(i, j), ROUNDED_ERROR);
			}
		}
	}
}

TEST_F(SvdOn7x5SparseMatrixTest, Multiplication_of_USV_should_be_input_matrix) {
	const double ROUNDED_ERROR = 1.0e-13;
	singular::Matrix< M, N > m2 =
		singular::Svd< M, N >::getU(this->usv)
		* singular::Svd< M, N >::getS(this->usv)
		* singular::Svd< M, N >::getV(this->usv).transpose();
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			EXPECT_NEAR(this->m(i, j), m2(i, j), ROUNDED_ERROR);
		}
	}
}

TEST_F(SvdOn7x5SparseMatrixTest, 5_positive_singular_values_should_be_produced) {
	const double ROUNDED_ERROR = 1.0e-14;
	const singular::DiagonalMatrix< M, N >& s =
		singular::Svd< M, N >::getS(this->usv);
	EXPECT_NEAR(14.0545840635182, s(0, 0), ROUNDED_ERROR * 10);
	EXPECT_NEAR(11.7483576435793, s(1, 1), ROUNDED_ERROR * 10);
	EXPECT_NEAR(11.5450620618378, s(2, 2), ROUNDED_ERROR * 10);
	EXPECT_NEAR(7.42393153794226, s(3, 3), ROUNDED_ERROR);
	EXPECT_NEAR(4.85402327850972, s(4, 4), ROUNDED_ERROR);
}

/** Fixture for SVD on a 2x2 matrix. */
class SvdOn2x2MatrixTest : public ::testing::Test {
protected:
	/** Number of rows in the input matrix. */
	static const int M = 2;

	/** Number of columns in the input matrix. */
	static const int N = 2;

	/** Input matrix. */
	singular::Matrix< M, N > m;

	/** Results of SVD. */
	singular::Svd< M, N >::USV usv;

	/** Builds the input matrix and performs SVD on it. */
	virtual void SetUp() {
		const double DATA[] = {
			0.0, 1.0,
			1.0, 0.0
		};
		this->m.fill(DATA);
		this->usv = singular::Svd< M, N >::decomposeUSV(this->m);
	}
};

TEST_F(SvdOn2x2MatrixTest, Left_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	singular::Matrix< M, M > eye =
		singular::Svd< M, N >::getU(this->usv)
		* singular::Svd< M, N >::getU(this->usv).transpose();
	EXPECT_NEAR(1.0, eye(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, eye(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, eye(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(1.0, eye(1, 1), ROUNDED_ERROR);
}

TEST_F(SvdOn2x2MatrixTest, Right_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	singular::Matrix< N, N > eye =
		singular::Svd< M, N >::getV(this->usv)
		* singular::Svd< M, N >::getV(this->usv).transpose();
	EXPECT_NEAR(1.0, eye(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, eye(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, eye(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(1.0, eye(1, 1), ROUNDED_ERROR);
}

TEST_F(SvdOn2x2MatrixTest, Multiplication_of_USV_should_be_input_matrix) {
	const double ROUNDED_ERROR = 1.0e-13;
	singular::Matrix< M, N > m2 =
		singular::Svd< M, N >::getU(this->usv)
		* singular::Svd< M, N >::getS(this->usv)
		* singular::Svd< M, N >::getV(this->usv).transpose();
	EXPECT_NEAR(this->m(0, 0), m2(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(this->m(0, 1), m2(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(this->m(1, 0), m2(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(this->m(1, 1), m2(1, 1), ROUNDED_ERROR);
}

TEST_F(SvdOn2x2MatrixTest, 2_positive_singular_values_should_be_produced) {
	const double ROUNDED_ERROR = 1.0e-14;
	const singular::DiagonalMatrix< M, N >& s =
		singular::Svd< M, N >::getS(this->usv);
	EXPECT_NEAR(1.0, s(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(1.0, s(0, 0), ROUNDED_ERROR);
}

/** Fixture for SVD on a 2x1 matrix. */
class SvdOn2x1MatrixTest : public ::testing::Test {
protected:
	/** Number of rows in the input matrix. */
	static const int M = 2;

	/** Number of columns in the input matrix. */
	static const int N = 1;

	/** Input matrix. */
	singular::Matrix< M, N > m;

	/** Results of SVD. */
	singular::Svd< M, N >::USV usv;

	/** Builds the input matrix and performs SVD on it. */
	virtual void SetUp() {
		const double DATA[] = {
			3.0,
			4.0
		};
		this->m.fill(DATA);
		this->usv = singular::Svd< M, N >::decomposeUSV(this->m);
	}
};

TEST_F(SvdOn2x1MatrixTest, Left_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	singular::Matrix< M, M > eye =
		singular::Svd< M, N >::getU(this->usv)
		* singular::Svd< M, N >::getU(this->usv).transpose();
	EXPECT_NEAR(1.0, eye(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, eye(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, eye(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(1.0, eye(1, 1), ROUNDED_ERROR);
}

TEST_F(SvdOn2x1MatrixTest, Right_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	singular::Matrix< N, N > eye =
		singular::Svd< M, N >::getV(this->usv)
		* singular::Svd< M, N >::getV(this->usv).transpose();
	EXPECT_NEAR(1.0, eye(0, 0), ROUNDED_ERROR);
}

TEST_F(SvdOn2x1MatrixTest, Multiplication_of_USV_should_be_input_matrix) {
	const double ROUNDED_ERROR = 1.0e-13;
	singular::Matrix< M, N > m2 =
		singular::Svd< M, N >::getU(this->usv)
		* singular::Svd< M, N >::getS(this->usv)
		* singular::Svd< M, N >::getV(this->usv).transpose();
	EXPECT_NEAR(this->m(0, 0), m2(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(this->m(1, 0), m2(1, 0), ROUNDED_ERROR);
}

TEST_F(SvdOn2x1MatrixTest, 1_positive_singular_value_should_be_produced) {
	const double ROUNDED_ERROR = 1.0e-14;
	const singular::DiagonalMatrix< M, N >& s =
		singular::Svd< M, N >::getS(this->usv);
	EXPECT_NEAR(5.0, s(0, 0), ROUNDED_ERROR);
}

/** Fixture for SVD on a 1x1 matrix. */
class SvdOn1x1MatrixTest : public ::testing::Test {
protected:
	/** Number of rows in the input matrix. */
	static const int M = 1;

	/** Number of columns in the input matrix. */
	static const int N = 1;

	/** Input matrix. */
	singular::Matrix< M, N > m;

	/** Results of SVD. */
	singular::Svd< M, N >::USV usv;

	/** Builds the input matrix and performs SVD on it. */
	virtual void SetUp() {
		const double DATA[] = {
			-1.0
		};
		this->m.fill(DATA);
		this->usv = singular::Svd< M, N >::decomposeUSV(this->m);
	}
};

TEST_F(SvdOn1x1MatrixTest, Left_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	singular::Matrix< M, M > eye =
		singular::Svd< M, N >::getU(this->usv)
		* singular::Svd< M, N >::getU(this->usv).transpose();
	EXPECT_NEAR(1.0, eye(0, 0), ROUNDED_ERROR);
}

TEST_F(SvdOn1x1MatrixTest, Right_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	singular::Matrix< N, N > eye =
		singular::Svd< M, N >::getV(this->usv)
		* singular::Svd< M, N >::getV(this->usv).transpose();
	EXPECT_NEAR(1.0, eye(0, 0), ROUNDED_ERROR);
}

TEST_F(SvdOn1x1MatrixTest, Multiplication_of_USV_should_be_input_matrix) {
	const double ROUNDED_ERROR = 1.0e-13;
	singular::Matrix< M, N > m2 =
		singular::Svd< M, N >::getU(this->usv)
		* singular::Svd< M, N >::getS(this->usv)
		* singular::Svd< M, N >::getV(this->usv).transpose();
	EXPECT_NEAR(this->m(0, 0), m2(0, 0), ROUNDED_ERROR);
}

TEST_F(SvdOn1x1MatrixTest, 1_positive_singular_value_should_be_produced) {
	const double ROUNDED_ERROR = 1.0e-14;
	const singular::DiagonalMatrix< M, N >& s =
		singular::Svd< M, N >::getS(this->usv);
	EXPECT_NEAR(1.0, s(0, 0), ROUNDED_ERROR);
}

/** Fixture for SVD on a 5x4 matrix filled with zeros. */
class SvdOn5x4ZerosTest : public ::testing::Test {
protected:
	/** Number of rows in the input matrix. */
	static const int M = 5;

	/** Number of columns in the input matrix. */
	static const int N = 4;

	/** Input matrix. */
	singular::Matrix< M, N > m;

	/** Results of SVD. */
	singular::Svd< M, N >::USV usv;

	/** Builds the input matrix and performs SVD on it. */
	virtual void SetUp() {
		const double DATA[] = {
			0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0
		};
		this->m.fill(DATA);
		this->usv = singular::Svd< M, N >::decomposeUSV(this->m);
	}
};

TEST_F(SvdOn5x4ZerosTest, Left_singular_matrix_should_be_identity_matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const singular::Matrix< M, M >& u = singular::Svd< M, N >::getU(this->usv);
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < M; ++j) {
			if (i == j) {
				EXPECT_NEAR(1.0, u(i, j), ROUNDED_ERROR);
			} else {
				EXPECT_NEAR(0.0, u(i, j), ROUNDED_ERROR);
			}
		}
	}
}

TEST_F(SvdOn5x4ZerosTest, Right_singular_matrix_should_be_identity_matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const singular::Matrix< N, N >& v = singular::Svd< M, N >::getV(this->usv);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (i == j) {
				EXPECT_NEAR(1.0, v(i, j), ROUNDED_ERROR);
			} else {
				EXPECT_NEAR(0.0, v(i, j), ROUNDED_ERROR);
			}
		}
	}
}

TEST_F(SvdOn5x4ZerosTest, Singular_values_should_be_zeros) {
	const double ROUNDED_ERROR = 1.0e-14;
	const singular::DiagonalMatrix< M, N >& s =
		singular::Svd< M, N >::getS(this->usv);
	for (int i = 0; i < N; ++i) {
		EXPECT_NEAR(0.0, s(i, i), ROUNDED_ERROR);
	}
}
