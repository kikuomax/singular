#include "singular/Vector.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <iterator>

TEST(VectorTest, Vector_can_be_associated_with_a_simple_array) {
	// [1.0, 2.0, 3.0, 4.0]
	double data[] = {
		1.0, 2.0, 3.0, 4.0
	};
	const size_t N = sizeof(data) / sizeof(double);
	singular::Vector< double > v(data, N, 1);
	// checks the size
	ASSERT_EQ(N, v.size());
	// checks the elements
	EXPECT_EQ(data[0], v[0]);
	EXPECT_EQ(data[1], v[1]);
	EXPECT_EQ(data[2], v[2]);
	EXPECT_EQ(data[3], v[3]);
	// checks the iterator
	ASSERT_EQ(N, std::distance(v.begin(), v.end()));
	EXPECT_TRUE(std::equal(v.begin(), v.end(), data));
}

TEST(VectorTest, Const_Vector_can_be_associated_with_a_simple_array) {
	// [1.0, 2.0, 3.0, 4.0]
	double data[] = {
		1.0, 2.0, 3.0, 4.0
	};
	const size_t N = sizeof(data) / sizeof(double);
	const singular::Vector< double > v(data, N, 1);
	// checks the size
	ASSERT_EQ(N, v.size());
	// checks the elements
	EXPECT_EQ(data[0], v[0]);
	EXPECT_EQ(data[1], v[1]);
	EXPECT_EQ(data[2], v[2]);
	EXPECT_EQ(data[3], v[3]);
	// checks the iterator
	ASSERT_EQ(N, std::distance(v.begin(), v.end()));
	EXPECT_TRUE(std::equal(v.begin(), v.end(), data));
}

TEST(VectorTest, Vector_can_be_associated_with_an_empty_array) {
	// empty data
	double data[1];  // data[0] is not allowed
	const size_t N = 0;
	singular::Vector< double > v(data, N, 1);
	// checks the size
	EXPECT_EQ(N, v.size());
	// checks the iterator
	EXPECT_EQ(N, std::distance(v.begin(), v.end()));
}

TEST(VectorTest, Vector_can_interleave_data_array) {
	// interleaved [1.5, 1.6, 1.7, 1.8]
	double block[] = {
		1.5, 0, 0,
		1.6, 0, 0,
		1.7, 0, 0,
		1.8, 0, 0
	};
	const double data[] = {
		1.5, 1.6, 1.7, 1.8
	};
	const size_t N = 4;
	singular::Vector< double > v(block, N, 3);
	// checks the size
	ASSERT_EQ(N, v.size());
	// checks the elements
	EXPECT_EQ(data[0], v[0]);
	EXPECT_EQ(data[1], v[1]);
	EXPECT_EQ(data[2], v[2]);
	EXPECT_EQ(data[3], v[3]);
	// checks the iterator
	ASSERT_EQ(N, std::distance(v.begin(), v.end()));
	EXPECT_TRUE(std::equal(v.begin(), v.end(), data));
}

TEST(VectorTest, Const_Vector_can_interleave_data_array) {
	// interleaved [1.5, 1.6, 1.7, 1.8]
	double block[] = {
		1.5, 0, 0,
		1.6, 0, 0,
		1.7, 0, 0,
		1.8, 0, 0
	};
	const double data[] = {
		1.5, 1.6, 1.7, 1.8
	};
	const size_t N = 4;
	const singular::Vector< double > v(block, N, 3);
	// checks the size
	ASSERT_EQ(N, v.size());
	// checks the elements
	EXPECT_EQ(data[0], v[0]);
	EXPECT_EQ(data[1], v[1]);
	EXPECT_EQ(data[2], v[2]);
	EXPECT_EQ(data[3], v[3]);
	// checks the iterator
	ASSERT_EQ(N, std::distance(v.begin(), v.end()));
	EXPECT_TRUE(std::equal(v.begin(), v.end(), data));
}

TEST(VectorTest, Vector_can_be_sliced) {
	double block[] = {
		1.0, 2.0, 3.0
	};
	singular::Vector< double > v(block, 3, 1);
	double data[] = {
		2.0, 3.0
	};
	const size_t N = 2;
	singular::Vector< double > v2 = v.slice(1);
	// checks the size
	ASSERT_EQ(N, v2.size());
	// checks the elements
	EXPECT_EQ(data[0], v2[0]);
	EXPECT_EQ(data[1], v2[1]);
	// checks the iterator
	ASSERT_EQ(N, std::distance(v2.begin(), v2.end()));
	EXPECT_TRUE(std::equal(v2.begin(), v2.end(), data));
}

TEST(VectorTest, Interleaved_Vector_can_be_sliced) {
	double block[] = {
		1.0, 0, 0,
		2.0, 0, 0,
		3.0, 0, 0,
		4.0, 0, 0
	};
	singular::Vector< double > v(block, 4, 3);
	const double data[] = {
		3.0, 4.0
	};
	const size_t N = 2;
	singular::Vector< double > v2 = v.slice(2);
	// checks the size
	ASSERT_EQ(N, v2.size());
	// checks the elements
	EXPECT_EQ(data[0], v2[0]);
	EXPECT_EQ(data[1], v2[1]);
	// checks the iterator
	ASSERT_EQ(N, std::distance(v2.begin(), v2.end()));
	EXPECT_TRUE(std::equal(v2.begin(), v2.end(), data));
}

TEST(VectorTest, Changes_on_Vector_should_be_reflected_to_source) {
	double block[] = {
		1.0, 2.0, 3.0
	};
	singular::Vector< double > v(block, 3, 1);
	v[0] = 0.0;
	v[1] = 1.5;
	v[2] = -1.0;
	EXPECT_EQ(v[0], 0.0);
	EXPECT_EQ(v[1], 1.5);
	EXPECT_EQ(v[2], -1.0);
	EXPECT_EQ(block[0], 0.0);
	EXPECT_EQ(block[1], 1.5);
	EXPECT_EQ(block[2], -1.0);
}

TEST(VectorTest, Changes_on_sliced_Vector_should_be_reflected_to_source) {
	double block[] = {
		1.5, 1.7, 2.6, 3.4
	};
	singular::Vector< double > v(block, 4, 1);
	singular::Vector< double > v2 = v.slice(2);
	v2[0] = -1.2;
	v2[1] = 0.9;
	EXPECT_EQ(v2[0], -1.2);
	EXPECT_EQ(v2[1], 0.9);
	EXPECT_EQ(v[0], 1.5);
	EXPECT_EQ(v[1], 1.7);
	EXPECT_EQ(v[2], -1.2);
	EXPECT_EQ(v[3], 0.9);
	EXPECT_EQ(block[0], 1.5);
	EXPECT_EQ(block[1], 1.7);
	EXPECT_EQ(block[2], -1.2);
	EXPECT_EQ(block[3], 0.9);
}
