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
	ASSERT_EQ(data[0], v[0]);
	ASSERT_EQ(data[1], v[1]);
	ASSERT_EQ(data[2], v[2]);
	ASSERT_EQ(data[3], v[3]);
	// checks the iterator
	ASSERT_EQ(N, std::distance(v.begin(), v.end()));
	ASSERT_TRUE(std::equal(v.begin(), v.end(), data));
}

TEST(VectorTest, Vector_can_be_associated_with_an_empty_array) {
	// empty data
	double data[0];
	const size_t N = 0;
	singular::Vector< double > v(data, N, 1);
	// checks the size
	ASSERT_EQ(N, v.size());
	// checks the iterator
	ASSERT_EQ(N, std::distance(v.begin(), v.end()));
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
	ASSERT_EQ(data[0], v[0]);
	ASSERT_EQ(data[1], v[1]);
	ASSERT_EQ(data[2], v[2]);
	ASSERT_EQ(data[3], v[3]);
	// checks the iterator
	ASSERT_EQ(N, std::distance(v.begin(), v.end()));
	ASSERT_TRUE(std::equal(v.begin(), v.end(), data));
}
