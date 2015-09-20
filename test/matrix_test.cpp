#include "singular/Matrix.h"

#include "gtest/gtest.h"

TEST(MatrixTest, 3x4_matrix_is_initally_fillwed_with_zeros) {
	const int M = 3;
	const int N = 4;
	singular::Matrix< M, N > m;
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			EXPECT_DOUBLE_EQ(0.0, m(i, j));
		}
	}
}

TEST(MatrixTest, 3x1_matrix_can_be_filled) {
	const int M = 3;
	const int N = 1;
	singular::Matrix< M, N > m;
	const double data[] = {
		1.0, 0.5, -2.0
	};
	singular::Matrix< M, N >& rM = m.fill(data);
	EXPECT_EQ(&m, &rM);
	EXPECT_DOUBLE_EQ(1.0, m(0, 0));
	EXPECT_DOUBLE_EQ(0.5, m(1, 0));
	EXPECT_DOUBLE_EQ(-2.0, m(2, 0));
}

TEST(MatrixTest, 1x4_matrix_can_be_filled) {
	const int M = 1;
	const int N = 4;
	singular::Matrix< M, N > m;
	const double data[] = {
		0.0, -5.0, 9.0, 1.0
	};
	singular::Matrix< M, N >& rM = m.fill(data);
	EXPECT_EQ(&m, &rM);
	EXPECT_DOUBLE_EQ(0.0, m(0, 0));
	EXPECT_DOUBLE_EQ(-5.0, m(0, 1));
	EXPECT_DOUBLE_EQ(9.0, m(0, 2));
	EXPECT_DOUBLE_EQ(1.0, m(0, 3));
}

TEST(MatrixTest, 3x4_matrix_can_be_filled) {
	const int M = 4;
	const int N = 3;
	singular::Matrix< M, N > m;
	const double data[] = {
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0,
		0.5, 1.5, 2.5
	};
	singular::Matrix< M, N >& rM = m.fill(data);
	EXPECT_EQ(&m, &rM);
	EXPECT_DOUBLE_EQ(1.0, m(0, 0));
	EXPECT_DOUBLE_EQ(2.0, m(0, 1));
	EXPECT_DOUBLE_EQ(3.0, m(0, 2));
	EXPECT_DOUBLE_EQ(4.0, m(1, 0));
	EXPECT_DOUBLE_EQ(5.0, m(1, 1));
	EXPECT_DOUBLE_EQ(6.0, m(1, 2));
	EXPECT_DOUBLE_EQ(7.0, m(2, 0));
	EXPECT_DOUBLE_EQ(8.0, m(2, 1));
	EXPECT_DOUBLE_EQ(9.0, m(2, 2));
	EXPECT_DOUBLE_EQ(0.5, m(3, 0));
	EXPECT_DOUBLE_EQ(1.5, m(3, 1));
	EXPECT_DOUBLE_EQ(2.5, m(3, 2));
}
