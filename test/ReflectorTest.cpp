#include "singular/Reflector.h"

#include "gtest/gtest.h"

TEST(ReflectorTest, Reflector_can_transform_a_4x1_matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const double DATA[] = {
		1.0, 2.0, 3.0, 4.0
	};
	singular::Matrix< 4, 1 > m = singular::Matrix< 4, 1 >::filledWith(DATA);
	singular::Reflector< 4 > h(singular::Vector< const double >(DATA, 4, 1));
	singular::Matrix< 4, 1 > m2 = h.applyFromLeftTo(m);
	EXPECT_NEAR(-5.477225575051661, m2(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(3, 0), ROUNDED_ERROR);
}

TEST(ReflectorTest, Reflector_can_transform_a_1x3_matrix) {
	const double ROUNDED_ERROR = 1.0e-14;
	const double DATA[] = {
		-1.0, -2.0, -3.0
	};
	singular::Matrix< 1, 3 > m = singular::Matrix< 1, 3 >::filledWith(DATA);
	singular::Reflector< 3 > h(singular::Vector< const double >(DATA, 3, 1));
	singular::Matrix< 1, 3 > m2 = h.applyFromRightTo(m);
	EXPECT_NEAR(3.741657386773941, m2(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(0, 2), ROUNDED_ERROR);
}

TEST(ReflectorTest, Reflector_can_transform_a_4x3_matrix_from_left) {
	const double ROUNDED_ERROR = 1.0e-14;
	const double DATA[] = {
		1.0,  2.0,  2.0,
		1.0,  0.5, -3.0,
		1.0, -2.0,  1.5,
		1.0,  3.0,  2.0
	};
	singular::Matrix< 4, 3 > m = singular::Matrix< 4, 3 >::filledWith(DATA);
	singular::Reflector< 4 > h(singular::Vector< const double >(DATA, 4, 3));
	singular::Matrix< 4, 3 > m2 = h.applyFromLeftTo(m);
	EXPECT_NEAR(-2.0, m2(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(-1.75, m2(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(-1.25, m2(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(-0.75, m2(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(-4.083333333333333, m2(1, 2), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(-3.25, m2(2, 1), ROUNDED_ERROR);
	EXPECT_NEAR(0.416666666666667, m2(2, 2), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(3, 0), ROUNDED_ERROR);
	EXPECT_NEAR(1.75, m2(3, 1), ROUNDED_ERROR);
	EXPECT_NEAR(0.916666666666667, m2(3, 2), ROUNDED_ERROR);
}

TEST(ReflectorTest, Reflector_can_transform_a_4x3_matrix_from_right) {
	const double ROUNDED_ERROR = 1.0e-14;
	const double DATA[] = {
		1.0,  2.0,  2.0,
		1.0,  0.5, -3.0,
		1.0, -2.0,  1.5,
		1.0,  3.0,  2.0
	};
	singular::Matrix< 4, 3 > m = singular::Matrix< 4, 3 >::filledWith(DATA);
	singular::Reflector< 3 > h(singular::Vector< const double >(DATA, 3, 1));
	singular::Matrix< 4, 3> m2 = h.applyFromRightTo(m);
	EXPECT_NEAR(-3.0, m2(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR(1.333333333333333, m2(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.666666666666667, m2(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(-2.833333333333333, m2(1, 2), ROUNDED_ERROR);
	EXPECT_NEAR(0, m2(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(-2.5, m2(2, 1), ROUNDED_ERROR);
	EXPECT_NEAR(1.0, m2(2, 2), ROUNDED_ERROR);
	EXPECT_NEAR(-3.666666666666667, m2(3, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.666666666666667, m2(3, 1), ROUNDED_ERROR);
	EXPECT_NEAR(-0.333333333333333, m2(3, 2), ROUNDED_ERROR);
}

TEST(ReflectorTest, 4x4_Reflector_can_transform_a_5x4_matrix_from_left) {
	const double ROUNDED_ERROR = 1.0e-14;
	const double DATA[] = {
		1.0, 2.0,  3.0,  4.0,
		2.0, 1.0,  2.0,  2.0,
		3.0, 1.0,  0.5, -3.0,
		4.0, 1.0, -2.0,  1.5,
		5.0, 1.0,  3.0,  2.0
	};
	singular::Matrix< 5, 4 > m = singular::Matrix< 5, 4 >::filledWith(DATA);
	singular::Reflector< 5 >
		h(singular::Vector< const double >(DATA + 5, 4, 4));
	singular::Matrix< 5, 4 > m2 = h.applyFromLeftTo(m);
	EXPECT_NEAR(1.0, m2(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(2.0, m2(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(3.0, m2(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR(4.0, m2(0, 3), ROUNDED_ERROR);
	EXPECT_NEAR(-7.0, m2(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(-2.0, m2(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(-1.75, m2(1, 2), ROUNDED_ERROR);
	EXPECT_NEAR(-1.25, m2(1, 3), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(2, 1), ROUNDED_ERROR);
	EXPECT_NEAR(-0.75, m2(2, 2), ROUNDED_ERROR);
	EXPECT_NEAR(-4.083333333333333, m2(2, 3), ROUNDED_ERROR);
	EXPECT_NEAR(1.0, m2(3, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(3, 1), ROUNDED_ERROR);
	EXPECT_NEAR(-3.25, m2(3, 2), ROUNDED_ERROR);
	EXPECT_NEAR(0.416666666666667, m2(3, 3), ROUNDED_ERROR);
	EXPECT_NEAR(2.0, m2(4, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(4, 1), ROUNDED_ERROR);
	EXPECT_NEAR(1.75, m2(4, 2), ROUNDED_ERROR);
	EXPECT_NEAR(0.916666666666667, m2(4, 3), ROUNDED_ERROR);
}

TEST(ReflectorTest, 3x3_Reflector_can_transform_a_5x4_matrix_from_right) {
	const double ROUNDED_ERROR = 1.0e-14;
	const double DATA[] = {
		1.0, 2.0,  3.0,  4.0,
		2.0, 1.0,  2.0,  2.0,
		3.0, 1.0,  0.5, -3.0,
		4.0, 1.0, -2.0,  1.5,
		5.0, 1.0,  3.0,  2.0
	};
	singular::Matrix< 5, 4 > m = singular::Matrix< 5, 4 >::filledWith(DATA);
	singular::Reflector< 4 >
		h(singular::Vector< const double >(DATA + 5, 3, 1));
	singular::Matrix< 5, 4 > m2 = h.applyFromRightTo(m);
	EXPECT_NEAR(1.0, m2(0, 0), ROUNDED_ERROR);
	EXPECT_NEAR(-5.333333333333333, m2(0, 1), ROUNDED_ERROR);
	EXPECT_NEAR(-0.666666666666667, m2(0, 2), ROUNDED_ERROR);
	EXPECT_NEAR(0.333333333333333, m2(0, 3), ROUNDED_ERROR);
	EXPECT_NEAR(2.0, m2(1, 0), ROUNDED_ERROR);
	EXPECT_NEAR(-3.0, m2(1, 1), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(1, 2), ROUNDED_ERROR);
	EXPECT_NEAR(0.0, m2(1, 3), ROUNDED_ERROR);
	EXPECT_NEAR(3.0, m2(2, 0), ROUNDED_ERROR);
	EXPECT_NEAR(1.333333333333333, m2(2, 1), ROUNDED_ERROR);
	EXPECT_NEAR(0.666666666666667, m2(2, 2), ROUNDED_ERROR);
	EXPECT_NEAR(-2.833333333333333, m2(2, 3), ROUNDED_ERROR);
	EXPECT_NEAR(4.0, m2(3, 0), ROUNDED_ERROR);
	EXPECT_NEAR(0, m2(3, 1), ROUNDED_ERROR);
	EXPECT_NEAR(-2.5, m2(3, 2), ROUNDED_ERROR);
	EXPECT_NEAR(1.0, m2(3, 3), ROUNDED_ERROR);
	EXPECT_NEAR(5.0, m2(4, 0), ROUNDED_ERROR);
	EXPECT_NEAR(-3.666666666666667, m2(4, 1), ROUNDED_ERROR);
	EXPECT_NEAR(0.666666666666667, m2(4, 2), ROUNDED_ERROR);
	EXPECT_NEAR(-0.333333333333333, m2(4, 3), ROUNDED_ERROR);
}
