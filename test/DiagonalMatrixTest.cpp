#include "singular/DiagonalMatrix.h"
#include "singular/Matrix.h"

#include "gtest/gtest.h"

TEST(DiagonalMatrixTest, 3x3_DiagonalMatrix_should_be_filled_with_zeros_by_default) {
	const int M = 3;
	const int N = 3;
	singular::DiagonalMatrix< M, N > m;
	EXPECT_EQ(0.0, m(0, 0));
	EXPECT_EQ(0.0, m(0, 1));
	EXPECT_EQ(0.0, m(0, 2));
	EXPECT_EQ(0.0, m(1, 0));
	EXPECT_EQ(0.0, m(1, 1));
	EXPECT_EQ(0.0, m(1, 2));
	EXPECT_EQ(0.0, m(2, 0));
	EXPECT_EQ(0.0, m(2, 1));
	EXPECT_EQ(0.0, m(2, 2));
}

TEST(DiagonalMatrixTest, 3x1_DiagonalMatrix_should_be_fillwed_with_zeros_by_default) {
	const int M = 3;
	const int N = 1;
	singular::DiagonalMatrix< M, N > m;
	EXPECT_EQ(0.0, m(0, 0));
	EXPECT_EQ(0.0, m(1, 0));
	EXPECT_EQ(0.0, m(2, 0));
}

TEST(DiagonalMatrixTest, 1x3_DiagonalMatrix_should_be_filled_with_zeros_by_default) {
	const int M = 1;
	const int N = 3;
	singular::DiagonalMatrix< M, N > m;
	EXPECT_EQ(0.0, m(0, 0));
	EXPECT_EQ(0.0, m(0, 1));
	EXPECT_EQ(0.0, m(0, 2));
}

TEST(DiagonalMatrixTest, 3x3_DiagonalMatrix_can_be_initialized_from_diagonal_elements) {
	const int M = 3;
	const int N = 3;
	const double DIAGONAL[] = {
		1.0, 2.0, 3.0
	};
	singular::DiagonalMatrix< M, N > m(DIAGONAL);
	EXPECT_EQ(1.0, m(0, 0));
	EXPECT_EQ(0.0, m(0, 1));
	EXPECT_EQ(0.0, m(0, 2));
	EXPECT_EQ(0.0, m(1, 0));
	EXPECT_EQ(2.0, m(1, 1));
	EXPECT_EQ(0.0, m(1, 2));
	EXPECT_EQ(0.0, m(2, 0));
	EXPECT_EQ(0.0, m(2, 1));
	EXPECT_EQ(3.0, m(2, 2));
}

TEST(DiagonalMatrixTest, 3x1_DiagonalMatrix_can_be_initialized_from_diagonal_element) {
	const int M = 3;
	const int N = 1;
	const double DIAGONAL[] = {
		-0.5
	};
	singular::DiagonalMatrix< M, N > m(DIAGONAL);
	EXPECT_EQ(-0.5, m(0, 0));
	EXPECT_EQ(0.0, m(1, 0));
	EXPECT_EQ(0.0, m(2, 0));
}

TEST(DiagonalMatrixTest, 1x3_DiagonalMatrix_can_be_initialized_from_diagonal_element) {
	const int M = 1;
	const int N = 3;
	const double DIAGONAL[] = {
		1.9
	};
	singular::DiagonalMatrix< M, N > m(DIAGONAL);
	EXPECT_EQ(1.9, m(0, 0));
	EXPECT_EQ(0.0, m(0, 1));
	EXPECT_EQ(0.0, m(0, 2));
}

TEST(DiagonalMatrixTest, Transposition_of_3x3_DiagonalMatrix_should_be_3x3_DiagonalMatrix) {
	const int M = 3;
	const int N = 3;
	const double DIAGONAL[] = {
		0.1, -2.5, 3.0
	};
	singular::DiagonalMatrix< M, N > m(DIAGONAL);
	singular::DiagonalMatrix< N, M > mT = m.transpose();
	EXPECT_EQ(0.1, mT(0, 0));
	EXPECT_EQ(0.0, mT(0, 1));
	EXPECT_EQ(0.0, mT(0, 2));
	EXPECT_EQ(0.0, mT(1, 0));
	EXPECT_EQ(-2.5, mT(1, 1));
	EXPECT_EQ(0.0, mT(1, 2));
	EXPECT_EQ(0.0, mT(2, 0));
	EXPECT_EQ(0.0, mT(2, 1));
	EXPECT_EQ(3.0, mT(2, 2));
}

TEST(DiagonalMatrixTest, Transposition_of_3x2_DiagonalMatrix_should_be_2x3_DiagonalMatrix) {
	const int M = 3;
	const int N = 2;
	const double DIAGONAL[] = {
		2.4, 1.5
	};
	singular::DiagonalMatrix< M, N > m(DIAGONAL);
	singular::DiagonalMatrix< N, M > mT = m.transpose();
	EXPECT_EQ(2.4, mT(0, 0));
	EXPECT_EQ(0.0, mT(0, 1));
	EXPECT_EQ(0.0, mT(0, 2));
	EXPECT_EQ(0.0, mT(1, 0));
	EXPECT_EQ(1.5, mT(1, 1));
	EXPECT_EQ(0.0, mT(1, 2));
}

TEST(DiagonalMatrixTest, Transposition_of_2x3_DiagonalMatrix_should_be_3x2_DiagonalMatrix) {
	const int M = 2;
	const int N = 3;
	const double DIAGONAL[] = {
		-1.2, -4.7
	};
	singular::DiagonalMatrix< M, N > m(DIAGONAL);
	singular::DiagonalMatrix< N, M > mT = m.transpose();
	EXPECT_EQ(-1.2, mT(0, 0));
	EXPECT_EQ(0.0, mT(0, 1));
	EXPECT_EQ(0.0, mT(1, 0));
	EXPECT_EQ(-4.7, mT(1, 1));
	EXPECT_EQ(0.0, mT(2, 0));
	EXPECT_EQ(0.0, mT(2, 1));
}

TEST(DiagonalMatrixTest, Product_of_3x4_DiagonalMatrix_and_4x3_Matrix_should_be_3x3_Matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M1 = 3;
	const int N1 = 4;
	const int M2 = 4;
	const int N2 = 3;
	const double DIAGONAL[] = {
		1.0, 2.0, 3.0
	};
	const double DATA[] = {
		3.5, 1.0, -2.4,
		1.5, 2.4, 1.0,
		-3.0, 1.9, 4.7,
		5.1, -3.2, 2.8
	};
	singular::DiagonalMatrix< M1, N1 > m1(DIAGONAL);
	singular::Matrix< M2, N2 > m2 =
		singular::Matrix< M2, N2 >::filledWith(DATA);
	singular::Matrix< M1, N2 > p = m1 * m2;
	EXPECT_NEAR(1.0 * 3.5, p(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(1.0 * 1.0, p(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(1.0 * -2.4, p(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR(2.0 * 1.5, p(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(2.0 * 2.4, p(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(2.0 * 1.0, p(1, 2), ROUNDED_ERROR);
	EXPECT_NEAR(3.0 * -3.0, p(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(3.0 * 1.9, p(2, 1), ROUNDED_ERROR);
	EXPECT_NEAR(3.0 * 4.7, p(2, 2), ROUNDED_ERROR);
}

TEST(DiagonalMatrixTest, Product_of_4x3_Matrix_and_3x4_DiagonalMatrix_should_be_4x4_Matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M1 = 4;
	const int N1 = 3;
	const int M2 = 3;
	const int N2 = 4;
	const double DATA[] = {
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0,
		10.0, 11.0, 12.0
	};
	const double DIAGONAL[] = {
		0.5, -2.1, 1.7
	};
	singular::Matrix< M1, N1 > m1 =
		singular::Matrix< M1, N1 >::filledWith(DATA);
	singular::DiagonalMatrix< M2, N2 > m2(DIAGONAL);
	singular::Matrix< M1, N2 > p = m1 * m2;
	EXPECT_NEAR(1.0 * 0.5, p(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(2.0 * -2.1, p(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(3.0 * 1.7, p(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, p(0, 3), ROUNDED_ERROR);
	EXPECT_NEAR(4.0 * 0.5, p(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(5.0 * -2.1, p(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(6.0 * 1.7, p(1, 2), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, p(1, 3), ROUNDED_ERROR);
	EXPECT_NEAR(7.0 * 0.5, p(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(8.0 * -2.1, p(2, 1), ROUNDED_ERROR);
	EXPECT_NEAR(9.0 * 1.7, p(2, 2), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, p(2, 3), ROUNDED_ERROR);
	EXPECT_NEAR(10.0 * 0.5, p(3, 0), ROUNDED_ERROR);
	EXPECT_NEAR(11.0 * -2.1, p(3, 1), ROUNDED_ERROR);
	EXPECT_NEAR(12.0 * 1.7, p(3, 2), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, p(3, 3), ROUNDED_ERROR);
}

TEST(DiagonalMatrixTest, Product_of_3x3_DiagonalMatrix_and_3x3_DiagonalMatrix_should_be_3x3_Matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M1 = 3;
	const int N1 = 3;
	const int M2 = 3;
	const int N2 = 3;
	const double DIAGONAL1[] = {
		1.0, 0.5, -3.2
	};
	const double DIAGONAL2[] = {
		-4.1, 2.9, 0.2
	};
	singular::DiagonalMatrix< M1, N1 > m1(DIAGONAL1);
	singular::DiagonalMatrix< M2, N2 > m2(DIAGONAL2);
	singular::Matrix< M1, N2 > p = m1 * m2;
	EXPECT_NEAR(1.0 * -4.1, p(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, p(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, p(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, p(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.5 * 2.9, p(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, p(1, 2), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, p(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, p(2, 1), ROUNDED_ERROR);
	EXPECT_NEAR(-3.2 * 0.2, p(2, 2), ROUNDED_ERROR);
}
