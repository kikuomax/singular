#include "singular/Matrix.h"

#include "gtest/gtest.h"

TEST(MatrixTest, 3x4_matrix_should_initally_be_fillwed_with_zeros) {
	const int M = 3;
	const int N = 4;
	singular::Matrix< M, N > m;
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			EXPECT_EQ(0.0, m(i, j));
		}
	}
}

TEST(MatrixTest, 3x3_identity_matrix_can_be_created) {
	const int M = 3;
	const int N = 3;
	singular::Matrix< M, N > eye = singular::Matrix< M, N >::identity();
	EXPECT_EQ(1.0, eye(0, 0));
	EXPECT_EQ(0.0, eye(0, 1));
	EXPECT_EQ(0.0, eye(0, 1));
	EXPECT_EQ(0.0, eye(1, 0));
	EXPECT_EQ(1.0, eye(1, 1));
	EXPECT_EQ(0.0, eye(1, 2));
	EXPECT_EQ(0.0, eye(2, 0));
	EXPECT_EQ(0.0, eye(2, 1));
	EXPECT_EQ(1.0, eye(2, 2));
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
	EXPECT_EQ(1.0, m(0, 0));
	EXPECT_EQ(0.5, m(1, 0));
	EXPECT_EQ(-2.0, m(2, 0));
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
	EXPECT_EQ(0.0, m(0, 0));
	EXPECT_EQ(-5.0, m(0, 1));
	EXPECT_EQ(9.0, m(0, 2));
	EXPECT_EQ(1.0, m(0, 3));
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
	EXPECT_EQ(1.0, m(0, 0));
	EXPECT_EQ(2.0, m(0, 1));
	EXPECT_EQ(3.0, m(0, 2));
	EXPECT_EQ(4.0, m(1, 0));
	EXPECT_EQ(5.0, m(1, 1));
	EXPECT_EQ(6.0, m(1, 2));
	EXPECT_EQ(7.0, m(2, 0));
	EXPECT_EQ(8.0, m(2, 1));
	EXPECT_EQ(9.0, m(2, 2));
	EXPECT_EQ(0.5, m(3, 0));
	EXPECT_EQ(1.5, m(3, 1));
	EXPECT_EQ(2.5, m(3, 2));
}

TEST(MatrixTest, Matrix_can_provide_row_vectors) {
	const int M = 4;
	const int N = 3;
	singular::Matrix< M, N > m;
	const double data[] = {
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0,
		0.5, 1.5, 2.5
	};
	m.fill(data);
	singular::Vector< const double > v = m.row(0);
	EXPECT_TRUE(std::equal(v.begin(), v.end(), data));
	v = m.row(1);
	EXPECT_TRUE(std::equal(v.begin(), v.end(), data + 3));
	v = m.row(2);
	EXPECT_TRUE(std::equal(v.begin(), v.end(), data + 6));
	v = m.row(3);
	EXPECT_TRUE(std::equal(v.begin(), v.end(), data + 9));
}

TEST(MatrixTest, Matrix_can_provide_column_vectors) {
	const int M = 4;
	const int N = 3;
	singular::Matrix< M, N > m;
	const double block[] = {
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0,
		0.5, 1.5, 2.5
	};
	m.fill(block);
	// column 0
	{
		const double data[] = { 1.0, 4.0, 7.0, 0.5 };
		singular::Vector< const double > v = m.column(0);
		EXPECT_TRUE(std::equal(v.begin(), v.end(), data));
	}
	// column 1
	{
		const double data[] = { 2.0, 5.0, 8.0, 1.5 };
		singular::Vector< const double > v = m.column(1);
		EXPECT_TRUE(std::equal(v.begin(), v.end(), data));
	}
	// column 2
	{
		const double data[] = { 3.0, 6.0, 9.0, 2.5 };
		singular::Vector< const double > v = m.column(2);
		EXPECT_TRUE(std::equal(v.begin(), v.end(), data));
	}
}
