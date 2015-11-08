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

TEST(MatrixTest, 3x3_matrix_can_be_cloned) {
	const int M = 3;
	const int N = 3;
	const double DATA[] = {
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0
	};
	singular::Matrix< M, N > m;
	m.fill(DATA);
	singular::Matrix< M, N > mC = m.clone();
	EXPECT_EQ(1.0, mC(0, 0));
	EXPECT_EQ(2.0, mC(0, 1));
	EXPECT_EQ(3.0, mC(0, 2));
	EXPECT_EQ(4.0, mC(1, 0));
	EXPECT_EQ(5.0, mC(1, 1));
	EXPECT_EQ(6.0, mC(1, 2));
	EXPECT_EQ(7.0, mC(2, 0));
	EXPECT_EQ(8.0, mC(2, 1));
	EXPECT_EQ(9.0, mC(2, 2));
}

TEST(MatrixTest, 3x4_matrix_can_be_cloned) {
	const int M = 3;
	const int N = 4;
	const double DATA[] = {
		1.0, 4.0, -7.0, 10.0,
		3.0, -6.0, 9.0, 12.0,
		-2.0, 5.0, 8.0, -11.0
	};
	singular::Matrix< M, N > m;
	m.fill(DATA);
	singular::Matrix< M, N > mC = m.clone();
	EXPECT_EQ(1.0, mC(0, 0));
	EXPECT_EQ(4.0, mC(0, 1));
	EXPECT_EQ(-7.0, mC(0, 2));
	EXPECT_EQ(10.0, mC(0, 3));
	EXPECT_EQ(3.0, mC(1, 0));
	EXPECT_EQ(-6.0, mC(1, 1));
	EXPECT_EQ(9.0, mC(1, 2));
	EXPECT_EQ(12.0, mC(1, 3));
	EXPECT_EQ(-2.0, mC(2, 0));
	EXPECT_EQ(5.0, mC(2, 1));
	EXPECT_EQ(8.0, mC(2, 2));
	EXPECT_EQ(-11.0, mC(2, 3));
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

TEST(MatrixTest, Product_of_3x3_matrix_and_3x3_matrix_should_be_3x3_matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M1 = 3;
	const int N1 = 3;
	const int M2 = 3;
	const int N2 = 3;
	const double DATA1[] = {
		1.0, -2.0, 3.0,
		3.0, 1.0, 2.0,
		4.0, 1.0, -5.0
	};
	const double DATA2[] = {
		2.0, 1.0, 4.0,
		5.0, 2.0, 3.0,
		6.0, -1.0, 3.0
	};
	singular::Matrix< M1, N1 > m1;
	m1.fill(DATA1);
	singular::Matrix< M2, N2 > m2;
	m2.fill(DATA2);
	singular::Matrix< M1, N2 > p = m1 * m2;
	EXPECT_NEAR(1.0*2.0 + -2.0*5.0 +  3.0* 6.0, p(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(1.0*1.0 + -2.0*2.0 +  3.0*-1.0, p(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(1.0*4.0 + -2.0*3.0 +  3.0* 3.0, p(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR(3.0*2.0 +  1.0*5.0 +  2.0* 6.0, p(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(3.0*1.0 +  1.0*2.0 +  2.0*-1.0, p(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(3.0*4.0 +  1.0*3.0 +  2.0* 3.0, p(1, 2), ROUNDED_ERROR);
	EXPECT_NEAR(4.0*2.0 +  1.0*5.0 + -5.0* 6.0, p(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(4.0*1.0 +  1.0*2.0 + -5.0*-1.0, p(2, 1), ROUNDED_ERROR);
	EXPECT_NEAR(4.0*4.0 +  1.0*3.0 + -5.0* 3.0, p(2, 2), ROUNDED_ERROR);
}

TEST(MatrixTest, Product_of_2x3_matrix_and_3x3_matrix_should_be_2x3_matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M1 = 2;
	const int N1 = 3;
	const int M2 = 3;
	const int N2 = 3;
	const double DATA1[] = {
		1.0, 2.0, 3.0,
		3.0, 1.0, 2.0
	};
	const double DATA2[] = {
		2.0, 5.0, -4.0,
		-1.0, 1.0, 1.0,
		4.0, 6.0, -3.0
	};
	singular::Matrix< M1, N1 > m1;
	m1.fill(DATA1);
	singular::Matrix< M2, N2 > m2;
	m2.fill(DATA2);
	singular::Matrix< M1, N2 > p = m1 * m2;
	EXPECT_NEAR(1.0* 2.0 + 2.0*-1.0 + 3.0* 4.0, p(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(1.0* 5.0 + 2.0* 1.0 + 3.0* 6.0, p(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(1.0*-4.0 + 2.0* 1.0 + 3.0*-3.0, p(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR(3.0* 2.0 + 1.0*-1.0 + 2.0* 4.0, p(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(3.0* 5.0 + 1.0* 1.0 + 2.0* 6.0, p(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(3.0*-4.0 + 1.0* 1.0 + 2.0*-3.0, p(1, 2), ROUNDED_ERROR);
}

TEST(MatrixTest, Product_of_4x3_matrix_and_3x3_matrix_should_be_4x3_matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M1 = 4;
	const int N1 = 3;
	const int M2 = 3;
	const int N2 = 3;
	const double DATA1[] = {
		1.0, 2.0, 3.0,
		4.0, 3.0, -1.0,
		2.0, -5.0, 4.0,
		-1.0, 2.0, 2.0
	};
	const double DATA2[] = {
		5.0, 2.0, 2.0,
		-1.0, 4.0, 3.0,
		2.0, 1.0, 6.0
	};
	singular::Matrix< M1, N1 > m1;
	m1.fill(DATA1);
	singular::Matrix< M2, N2 > m2;
	m2.fill(DATA2);
	singular::Matrix< M1, N2 > p = m1 * m2;
	EXPECT_NEAR( 1.0*5.0 +  2.0*-1.0 +  3.0*2.0, p(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR( 1.0*2.0 +  2.0* 4.0 +  3.0*1.0, p(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR( 1.0*2.0 +  2.0* 3.0 +  3.0*6.0, p(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR( 4.0*5.0 +  3.0*-1.0 + -1.0*2.0, p(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR( 4.0*2.0 +  3.0* 4.0 + -1.0*1.0, p(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR( 4.0*2.0 +  3.0* 3.0 + -1.0*6.0, p(1, 2), ROUNDED_ERROR);
	EXPECT_NEAR( 2.0*5.0 + -5.0*-1.0 +  4.0*2.0, p(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR( 2.0*2.0 + -5.0* 4.0 +  4.0*1.0, p(2, 1), ROUNDED_ERROR);
	EXPECT_NEAR( 2.0*2.0 + -5.0* 3.0 +  4.0*6.0, p(2, 2), ROUNDED_ERROR);
	EXPECT_NEAR(-1.0*5.0 +  2.0*-1.0 +  2.0*2.0, p(3, 0), ROUNDED_ERROR);
	EXPECT_NEAR(-1.0*2.0 +  2.0* 4.0 +  2.0*1.0, p(3, 1), ROUNDED_ERROR);
	EXPECT_NEAR(-1.0*2.0 +  2.0* 3.0 +  2.0*6.0, p(3, 2), ROUNDED_ERROR);
}

TEST(MatrixTest, Product_of_3x2_matrix_and_2x3_matrix_should_be_3x3_matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M1 = 3;
	const int N1 = 2;
	const int M2 = 2;
	const int N2 = 3;
	const double DATA1[] = {
		1.0, 2.0,
		3.0, 4.0,
		5.0, 6.0
	};
	const double DATA2[] = {
		-1.0, 5.0, 4.0,
		2.0, 1.0, -2.0
	};
	singular::Matrix< M1, N1 > m1;
	m1.fill(DATA1);
	singular::Matrix< M2, N2 > m2;
	m2.fill(DATA2);
	singular::Matrix< M1, N2 > p = m1 * m2;
	EXPECT_NEAR(1.0*-1.0 + 2.0* 2.0, p(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(1.0* 5.0 + 2.0* 1.0, p(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(1.0* 4.0 + 2.0*-2.0, p(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR(3.0*-1.0 + 4.0* 2.0, p(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(3.0* 5.0 + 4.0* 1.0, p(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(3.0* 4.0 + 4.0*-2.0, p(1, 2), ROUNDED_ERROR);
	EXPECT_NEAR(5.0*-1.0 + 6.0* 2.0, p(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(5.0* 5.0 + 6.0* 1.0, p(2, 1), ROUNDED_ERROR);
	EXPECT_NEAR(5.0* 4.0 + 6.0*-2.0, p(2, 2), ROUNDED_ERROR);
}

TEST(MatrixTest, Product_of_3x4_matrix_and_4x3_matrix_should_be_3x3_matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M1 = 3;
	const int N1 = 4;
	const int M2 = 4;
	const int N2 = 3;
	const double DATA1[] = {
		1.0, -2.0, 6.0, 3.0,
		4.0, 5.0, -1.0, 2.0,
		2.0, 3.0, -4.0, 1.0
	};
	const double DATA2[] = {
		2.0, -5.0, 3.0,
		1.0, 1.0, 2.0,
		4.0, 3.0, 6.0,
		2.0, 1.0, -2.0
	};
	singular::Matrix< M1, N1 > m1;
	m1.fill(DATA1);
	singular::Matrix< M2, N2 > m2;
	m2.fill(DATA2);
	singular::Matrix< M1, N2 > p = m1 * m2;
	EXPECT_NEAR(1.0* 2.0 + -2.0*1.0 +  6.0*4.0 + 3.0* 2.0, p(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(1.0*-5.0 + -2.0*1.0 +  6.0*3.0 + 3.0* 1.0, p(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(1.0* 3.0 + -2.0*2.0 +  6.0*6.0 + 3.0*-2.0, p(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR(4.0* 2.0 +  5.0*1.0 + -1.0*4.0 + 2.0* 2.0, p(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(4.0*-5.0 +  5.0*1.0 + -1.0*3.0 + 2.0* 1.0, p(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(4.0* 3.0 +  5.0*2.0 + -1.0*6.0 + 2.0*-2.0, p(1, 2), ROUNDED_ERROR);
	EXPECT_NEAR(2.0* 2.0 +  3.0*1.0 + -4.0*4.0 + 1.0* 2.0, p(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(2.0*-5.0 +  3.0*1.0 + -4.0*3.0 + 1.0* 1.0, p(2, 1), ROUNDED_ERROR);
	EXPECT_NEAR(2.0* 3.0 +  3.0*2.0 + -4.0*6.0 + 1.0*-2.0, p(2, 2), ROUNDED_ERROR);
}

TEST(MatrixTest, Product_of_3x3_matrix_and_3x2_matrix_should_be_3x2_matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M1 = 3;
	const int N1 = 3;
	const int M2 = 3;
	const int N2 = 2;
	const double DATA1[] = {
		6.0, 2.0, -3.0,
		1.0, 5.0, 4.0,
		4.0, 3.0, -1.0
	};
	const double DATA2[] = {
		3.0, 2.0,
		1.0, -2.0,
		4.0, 5.0
	};
	singular::Matrix< M1, N1 > m1;
	m1.fill(DATA1);
	singular::Matrix< M2, N2 > m2;
	m2.fill(DATA2);
	singular::Matrix< M1, N2 > p = m1 * m2;
	EXPECT_NEAR(6.0*3.0 + 2.0* 1.0 + -3.0*4.0, p(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(6.0*2.0 + 2.0*-2.0 + -3.0*5.0, p(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(1.0*3.0 + 5.0* 1.0 +  4.0*4.0, p(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(1.0*2.0 + 5.0*-2.0 +  4.0*5.0, p(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(4.0*3.0 + 3.0* 1.0 + -1.0*4.0, p(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(4.0*2.0 + 3.0*-2.0 + -1.0*5.0, p(2, 1), ROUNDED_ERROR);
}

TEST(MatrixTest, Product_of_3x3_matrix_and_3x4_matrix_should_be_3x4_matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M1 = 3;
	const int N1 = 3;
	const int M2 = 3;
	const int N2 = 4;
	const double DATA1[] = {
		3.0, -2.0, 1.0,
		2.0, -1.0, 4.0,
		-1.0, -3.0, 5.0
	};
	const double DATA2[] = {
		-4.0, 1.0, 5.0, 3.0,
		1.0, 2.0, 1.0, 2.0,
		-2.0, 3.0, 6.0, 1.0
	};
	singular::Matrix< M1, N1 > m1;
	m1.fill(DATA1);
	singular::Matrix< M2, N2 > m2;
	m2.fill(DATA2);
	singular::Matrix< M1, N2 > p = m1 * m2;
	EXPECT_NEAR( 3.0*-4.0 + -2.0*1.0 + 1.0*-2.0, p(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR( 3.0* 1.0 + -2.0*2.0 + 1.0* 3.0, p(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR( 3.0* 5.0 + -2.0*1.0 + 1.0* 6.0, p(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR( 3.0* 3.0 + -2.0*2.0 + 1.0* 1.0, p(0, 3), ROUNDED_ERROR);
	EXPECT_NEAR( 2.0*-4.0 + -1.0*1.0 + 4.0*-2.0, p(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR( 2.0* 1.0 + -1.0*2.0 + 4.0* 3.0, p(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR( 2.0* 5.0 + -1.0*1.0 + 4.0* 6.0, p(1, 2), ROUNDED_ERROR);
	EXPECT_NEAR( 2.0* 3.0 + -1.0*2.0 + 4.0* 1.0, p(1, 3), ROUNDED_ERROR);
	EXPECT_NEAR(-1.0*-4.0 + -3.0*1.0 + 5.0*-2.0, p(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(-1.0* 1.0 + -3.0*2.0 + 5.0* 3.0, p(2, 1), ROUNDED_ERROR);
	EXPECT_NEAR(-1.0* 5.0 + -3.0*1.0 + 5.0* 6.0, p(2, 2), ROUNDED_ERROR);
	EXPECT_NEAR(-1.0* 3.0 + -3.0*2.0 + 5.0* 1.0, p(2, 3), ROUNDED_ERROR);
}

TEST(MatrixTest, Transposition_of_3x3_matrix_should_be_3x3_matrix) {
	const int M = 3;
	const int N = 3;
	const double DATA[] = {
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0
	};
	singular::Matrix< M, N > m;
	m.fill(DATA);
	singular::Matrix< N, M > mT = m.transpose();
	EXPECT_EQ(1.0, mT(0, 0));
	EXPECT_EQ(4.0, mT(0, 1));
	EXPECT_EQ(7.0, mT(0, 2));
	EXPECT_EQ(2.0, mT(1, 0));
	EXPECT_EQ(5.0, mT(1, 1));
	EXPECT_EQ(8.0, mT(1, 2));
	EXPECT_EQ(3.0, mT(2, 0));
	EXPECT_EQ(6.0, mT(2, 1));
	EXPECT_EQ(9.0, mT(2, 2));
}

TEST(MatrixTest, Transposition_of_3x2_matrix_should_be_2x3_matrix) {
	const int M = 3;
	const int N = 2;
	const double DATA[] = {
		1.0, 4.0,
		2.0, 5.0,
		3.0, 6.0
	};
	singular::Matrix< M, N > m;
	m.fill(DATA);
	singular::Matrix< N, M > mT = m.transpose();
	EXPECT_EQ(1.0, mT(0, 0));
	EXPECT_EQ(2.0, mT(0, 1));
	EXPECT_EQ(3.0, mT(0, 2));
	EXPECT_EQ(4.0, mT(1, 0));
	EXPECT_EQ(5.0, mT(1, 1));
	EXPECT_EQ(6.0, mT(1, 2));
}

TEST(MatrixTest, Transposition_of_2x3_matrix_should_be_3x2_matrix) {
	const int M = 2;
	const int N = 3;
	const double DATA[] = {
		1.0, -2.0, 3.0,
		-4.0, 5.0, -6.0
	};
	singular::Matrix< M, N > m;
	m.fill(DATA);
	singular::Matrix< N, M > mT = m.transpose();
	EXPECT_EQ(1.0, mT(0, 0));
	EXPECT_EQ(-4.0, mT(0, 1));
	EXPECT_EQ(-2.0, mT(1, 0));
	EXPECT_EQ(5.0, mT(1, 1));
	EXPECT_EQ(3.0, mT(2, 0));
	EXPECT_EQ(-6.0, mT(2, 1));
}

TEST(MatrixTest, Rows_of_3x3_matrix_can_be_shuffled) {
	const int M = 3;
	const int N = 3;
	const double DATA[] = {
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0
	};
	const int ORDER[] = { 2, 1, 0 };
	singular::Matrix< M, N > m;
	m.fill(DATA);
	singular::Matrix< M, N > mS = m.shuffleRows(ORDER);
	EXPECT_EQ(7.0, mS(0, 0));
	EXPECT_EQ(8.0, mS(0, 1));
	EXPECT_EQ(9.0, mS(0, 2));
	EXPECT_EQ(4.0, mS(1, 0));
	EXPECT_EQ(5.0, mS(1, 1));
	EXPECT_EQ(6.0, mS(1, 2));
	EXPECT_EQ(1.0, mS(2, 0));
	EXPECT_EQ(2.0, mS(2, 1));
	EXPECT_EQ(3.0, mS(2, 2));
}

TEST(MatrixTest, Rows_of_4x3_matrix_can_be_shuffled) {
	const int M = 4;
	const int N = 3;
	const double DATA[] = {
		1.0, 5.0, 9.0,
		2.0, 6.0, 10.0,
		3.0, 7.0, 11.0,
		4.0, 8.0, 12.0
	};
	const int ORDER[] = { 0, 3, 1, 2 };
	singular::Matrix< M, N > m;
	m.fill(DATA);
	singular::Matrix< M, N > mS = m.shuffleRows(ORDER);
	EXPECT_EQ(1.0, mS(0, 0));
	EXPECT_EQ(5.0, mS(0, 1));
	EXPECT_EQ(9.0, mS(0, 2));
	EXPECT_EQ(4.0, mS(1, 0));
	EXPECT_EQ(8.0, mS(1, 1));
	EXPECT_EQ(12.0, mS(1, 2));
	EXPECT_EQ(2.0, mS(2, 0));
	EXPECT_EQ(6.0, mS(2, 1));
	EXPECT_EQ(10.0, mS(2, 2));
	EXPECT_EQ(3.0, mS(3, 0));
	EXPECT_EQ(7.0, mS(3, 1));
	EXPECT_EQ(11.0, mS(3, 2));
}

TEST(MatrixTest, Rows_of_3x4_matrix_can_be_shuffled) {
	const int M = 3;
	const int N = 4;
	const double DATA[] = {
		1.0, 3.0, 5.0, 7.0,
		8.0, 6.0, 4.0, 2.0,
		9.0, 10.0, 11.0, 12.0
	};
	const int ORDER[] = { 1, 2, 0 };
	singular::Matrix< M, N > m;
	m.fill(DATA);
	singular::Matrix< M, N > mS = m.shuffleRows(ORDER);
	EXPECT_EQ(8.0, mS(0, 0));
	EXPECT_EQ(6.0, mS(0, 1));
	EXPECT_EQ(4.0, mS(0, 2));
	EXPECT_EQ(2.0, mS(0, 3));
	EXPECT_EQ(9.0, mS(1, 0));
	EXPECT_EQ(10.0, mS(1, 1));
	EXPECT_EQ(11.0, mS(1, 2));
	EXPECT_EQ(12.0, mS(1, 3));
	EXPECT_EQ(1.0, mS(2, 0));
	EXPECT_EQ(3.0, mS(2, 1));
	EXPECT_EQ(5.0, mS(2, 2));
	EXPECT_EQ(7.0, mS(2, 3));
}

TEST(MatrixTest, Columns_of_3x3_matrix_can_be_shuffled) {
	const int M = 3;
	const int N = 3;
	const double DATA[] = {
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0
	};
	const int ORDER[] = { 2, 1, 0 };
	singular::Matrix< M, N > m;
	m.fill(DATA);
	singular::Matrix< M, N > mS = m.shuffleColumns(ORDER);
	EXPECT_EQ(3.0, mS(0, 0));
	EXPECT_EQ(6.0, mS(1, 0));
	EXPECT_EQ(9.0, mS(2, 0));
	EXPECT_EQ(2.0, mS(0, 1));
	EXPECT_EQ(5.0, mS(1, 1));
	EXPECT_EQ(8.0, mS(2, 1));
	EXPECT_EQ(1.0, mS(0, 2));
	EXPECT_EQ(4.0, mS(1, 2));
	EXPECT_EQ(7.0, mS(2, 2));
}

TEST(MatrixTest, Columns_of_4x3_matrix_can_be_shuffled) {
	const int M = 4;
	const int N = 3;
	const double DATA[] = {
		1.0, 5.0, 9.0,
		2.0, 6.0, 10.0,
		3.0, 7.0, 11.0,
		4.0, 8.0, 12.0
	};
	const int ORDER[] = { 1, 2, 0 };
	singular::Matrix< M, N > m;
	m.fill(DATA);
	singular::Matrix< M, N > mS = m.shuffleColumns(ORDER);
	EXPECT_EQ(5.0, mS(0, 0));
	EXPECT_EQ(6.0, mS(1, 0));
	EXPECT_EQ(7.0, mS(2, 0));
	EXPECT_EQ(8.0, mS(3, 0));
	EXPECT_EQ(9.0, mS(0, 1));
	EXPECT_EQ(10.0, mS(1, 1));
	EXPECT_EQ(11.0, mS(2, 1));
	EXPECT_EQ(12.0, mS(3, 1));
	EXPECT_EQ(1.0, mS(0, 2));
	EXPECT_EQ(2.0, mS(1, 2));
	EXPECT_EQ(3.0, mS(2, 2));
	EXPECT_EQ(4.0, mS(3, 2));
}

TEST(MatrixTest, Columns_of_3x4_matrix_can_be_shuffled) {
	const int M = 3;
	const int N = 4;
	const double DATA[] = {
		1.0, 3.0, 5.0, 7.0,
		8.0, 6.0, 4.0, 2.0,
		9.0, 10.0, 11.0, 12.0
	};
	const int ORDER[] = { 0, 3, 1, 2 };
	singular::Matrix< M, N > m;
	m.fill(DATA);
	singular::Matrix< M, N > mS = m.shuffleColumns(ORDER);
	EXPECT_EQ(1.0, mS(0, 0));
	EXPECT_EQ(8.0, mS(1, 0));
	EXPECT_EQ(9.0, mS(2, 0));
	EXPECT_EQ(7.0, mS(0, 1));
	EXPECT_EQ(2.0, mS(1, 1));
	EXPECT_EQ(12.0, mS(2, 1));
	EXPECT_EQ(3.0, mS(0, 2));
	EXPECT_EQ(6.0, mS(1, 2));
	EXPECT_EQ(10.0, mS(2, 2));
	EXPECT_EQ(5.0, mS(0, 3));
	EXPECT_EQ(4.0, mS(1, 3));
	EXPECT_EQ(11.0, mS(2, 3));
}
