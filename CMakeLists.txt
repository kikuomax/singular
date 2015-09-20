cmake_minimum_required (VERSION 2.8.4)

project (singular)

set (PROJECT_MAJOR 0)
set (PROJECT_MINOR 1)
set (PROJECT_PATCH 0-SNAPSHOT)
set (PROJECT_VERSION "${PROJECT_MAJOR}.${PROJECT_MINOR}.${PROJECT_PATCH}")

# enables C++11 on GCC like compilers
if (CMAKE_COMPILER_IS_GNUCXX OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif ()

# finds Google Test
find_package (GTest)

# turn on if you want to generate tests
# turned on by default as far as GTest is found
option (GENERATE_TESTS "Generate Google Test tests" ${GTEST_FOUND})

# finds doxygen
find_package (Doxygen)

# turn on if you want to generate documentation
# turned on by default as far as Doxygen is found
option (GENERATE_DOCUMENTATION
	"Generate Doxygen documentation" ${DOXYGEN_FOUND})

include_directories (
	src
	${GTEST_INCLUDE_DIRS})

add_executable (svd-test test/test.cpp)

# generates tests if necessary
if (GENERATE_TESTS)
	enable_testing ()

	add_executable (singular-test
		test/vector_test.cpp
		test/matrix_test.cpp
		test/reflector_test.cpp)

	target_link_libraries (singular-test ${GTEST_BOTH_LIBRARIES})

	add_test (NAME singular-test COMMAND singular-test)
endif ()

# generates documentation if necessary
if (GENERATE_DOCUMENTATION)
	# configures Doxygen
	configure_file (Doxyfile.in Doxyfile @ONLY)

	add_custom_target (docs ALL
		${DOXYGEN_EXECUTABLE} Doxyfile
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
		COMMENT "Generating Doxygen documentation")
endif ()
