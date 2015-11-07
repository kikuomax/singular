#include "singular/svd.h"

#include "gtest/gtest.h"

TEST(SVDTest, 5x5_left_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M = 5;
	const int N = 4;
	const double DATA[] = {
		1.0, 2.0, 3.0, 4.0,
		5.0, 6.0, 7.0, 8.0,
		4.0, 8.0, 3.0, 5.0,
		6.0, 7.0, 2.0, 1.0,
		9.0, 1.0, 3.0, 6.0
	};
	singular::Matrix< M, N > m;
	m.fill(DATA);
	singular::Svd< M, N >::USV usv = singular::Svd< M, N >::decomposeUSV(m);
	singular::Matrix< M, M > eye =
		singular::Svd< M, N >::getU(usv)
		* singular::Svd< M, N >::getU(usv).transpose();
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

TEST(SVDTest, 4x4_right_singular_vectors_should_be_orthonormal) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M = 5;
	const int N = 4;
	const double DATA[] = {
		1.0, 2.0, 3.0, 4.0,
		5.0, 6.0, 7.0, 8.0,
		4.0, 8.0, 3.0, 5.0,
		6.0, 7.0, 2.0, 1.0,
		9.0, 1.0, 3.0, 6.0
	};
	singular::Matrix< M, N > m;
	m.fill(DATA);
	singular::Svd< M, N >::USV usv = singular::Svd< M, N >::decomposeUSV(m);
	singular::Matrix< N, N > eye =
		singular::Svd< M, N >::getV(usv)
		* singular::Svd< M, N >::getV(usv).transpose();
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

TEST(SVDTest, Multiplication_of_USV_should_be_original_5x4_matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M = 5;
	const int N = 4;
	const double DATA[] = {
		1.0, 2.0, 3.0, 4.0,
		5.0, 6.0, 7.0, 8.0,
		4.0, 8.0, 3.0, 5.0,
		6.0, 7.0, 2.0, 1.0,
		9.0, 1.0, 3.0, 6.0
	};
	singular::Matrix< M, N > m;
	m.fill(DATA);
	singular::Svd< M, N >::USV usv = singular::Svd< M, N >::decomposeUSV(m);
	singular::Matrix< M, N > m2 =
		singular::Svd< M, N >::getU(usv)
		* singular::Svd< M, N >::getS(usv)
		* singular::Svd< M, N >::getV(usv).transpose();
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			EXPECT_NEAR(m(i, j), m2(i, j), ROUNDED_ERROR);
		}
	}
}
