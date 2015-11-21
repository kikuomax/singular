#include "singular/Rotator.h"

#include "gtest/gtest.h"

TEST(RotatorTest, Rotator_should_be_a_2x2_matrix) {
	const double ROUNDED_ERROR = 1.0e-15;
	// rotator for [1, 1]
	{
		singular::Rotator r(1, 1);
		EXPECT_NEAR(0.707106781186548, r(0, 0), ROUNDED_ERROR);
		EXPECT_NEAR(-0.707106781186548, r(0, 1), ROUNDED_ERROR);
		EXPECT_NEAR(0.707106781186548, r(1, 0), ROUNDED_ERROR);
		EXPECT_NEAR(0.707106781186548, r(1, 1), ROUNDED_ERROR);
	}
	// rotator for [3, -2]
	{
		singular::Rotator r(3, -2);
		EXPECT_NEAR(0.832050294337844, r(0, 0), ROUNDED_ERROR);
		EXPECT_NEAR(0.554700196225229, r(0, 1), ROUNDED_ERROR);
		EXPECT_NEAR(-0.554700196225229, r(1, 0), ROUNDED_ERROR);
		EXPECT_NEAR(0.832050294337844, r(1, 1), ROUNDED_ERROR);
	}
	// rotator for [1, 0]
	{
		singular::Rotator r(1, 0);
		EXPECT_NEAR(1.0, r(0, 0), ROUNDED_ERROR);
		EXPECT_NEAR(0.0, r(0, 1), ROUNDED_ERROR);
		EXPECT_NEAR(0.0, r(1, 0), ROUNDED_ERROR);
		EXPECT_NEAR(1.0, r(1, 1), ROUNDED_ERROR);
	}
	// rotator for [0, 1]
	{
		singular::Rotator r(0, 1);
		EXPECT_NEAR(0.0, r(0, 0), ROUNDED_ERROR);
		EXPECT_NEAR(-1.0, r(0, 1), ROUNDED_ERROR);
		EXPECT_NEAR(1.0, r(1, 0), ROUNDED_ERROR);
		EXPECT_NEAR(0.0, r(1, 1), ROUNDED_ERROR);
	}
}

TEST(RotatorTest, Rotator_can_transform_2x1_matrix_from_left) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M = 2;
	const int N = 1;
	const double DATA[] = { 1, 2 };
	singular::Matrix< M, N > m = singular::Matrix< M, N >::filledWith(DATA);
	singular::Rotator r(1, 2);
	singular::Matrix< M, N > m2 = r.applyFromLeftTo(m, 0);
	EXPECT_NEAR(2.236067977499790, m2(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(1, 0), ROUNDED_ERROR);
}

TEST(RotatorTest, Rotator_can_transform_1x2_matrix_from_right) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M = 1;
	const int N = 2;
	const double DATA[] = { 1, 2 };
	singular::Matrix< M, N > m = singular::Matrix< M, N >::filledWith(DATA);
	singular::Rotator r(1, 2);
	singular::Matrix< M, N > m2 = r.applyFromRightTo(m, 0);
	EXPECT_NEAR(2.236067977499790, m2(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(0, 1), ROUNDED_ERROR);
}

TEST(RotatorTest, Rotator_can_transform_4x3_matrix_from_left) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M = 4;
	const int N = 3;
	const double DATA[] = {
		1, 3, 8,
		2, 6, 5,
		4, 2, 7,
		8, 9, 1
	};
	singular::Matrix< M, N > m = singular::Matrix< M, N >::filledWith(DATA);
	singular::Rotator r(6, 2);
	singular::Matrix< M, N > m2 = r.applyFromLeftTo(m, 1);
	EXPECT_NEAR(1, m2(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(3, m2(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(8, m2(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR(3.162277660168379, m2(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(6.324555320336759, m2(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(6.957010852370435, m2(1, 2), ROUNDED_ERROR);
	EXPECT_NEAR(3.162277660168379, m2(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0, m2(2, 1), ROUNDED_ERROR);
	EXPECT_NEAR(5.059644256269407, m2(2, 2), ROUNDED_ERROR);
	EXPECT_NEAR(8, m2(3, 0), ROUNDED_ERROR);
	EXPECT_NEAR(9, m2(3, 1), ROUNDED_ERROR);
	EXPECT_NEAR(1, m2(3, 2), ROUNDED_ERROR);
}

TEST(RotatorTest, Rotator_can_transform_4x3_matrix_from_right) {
	const double ROUNDED_ERROR = 1.0e-14;
	const int M = 4;
	const int N = 3;
	const double DATA[] = {
		1, 3, 8,
		2, 6, 5,
		4, 2, 7,
		8, 9, 1
	};
	singular::Matrix< M, N > m = singular::Matrix< M, N >::filledWith(DATA);
	singular::Rotator r(6, 5);
	singular::Matrix< M, N > m2 = r.applyFromRightTo(m, 1);
	EXPECT_NEAR(1, m2(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(7.426139036107967, m2(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(4.225217037785567, m2(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR(2, m2(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(7.810249675906654, m2(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(0, m2(1, 2), ROUNDED_ERROR);
	EXPECT_NEAR(4, m2(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(6.017733356846111, m2(2, 1), ROUNDED_ERROR);
	EXPECT_NEAR(4.097180157852671, m2(2, 2), ROUNDED_ERROR);
	EXPECT_NEAR(8, m2(3, 0), ROUNDED_ERROR);
	EXPECT_NEAR(7.554175916040862, m2(3, 1), ROUNDED_ERROR);
	EXPECT_NEAR(-4.993438317382943, m2(3, 2), ROUNDED_ERROR);
}
