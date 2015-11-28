cmake_minimum_required (VERSION 2.8.4)

project (singular)

set (PROJECT_MAJOR 0)
set (PROJECT_MINOR 1)
set (PROJECT_PATCH 1)
set (PROJECT_VERSION "${PROJECT_MAJOR}.${PROJECT_MINOR}.${PROJECT_PATCH}")

# turn off if you want to dynamically link the MSVC runtime
option (MSVC_LINK_STATIC_RUNTIME "Statically link the MSVC runtime" on)

# replaces /MD --> /MT if the static MSVC runtime is needed
if (MSVC AND MSVC_LINK_STATIC_RUNTIME)
	foreach (flag_var
		CMAKE_CXX_FLAGS
		CMAKE_CXX_FLAGS_DEBUG
		CMAKE_CXX_FLAGS_RELEASE
		CMAKE_CXX_FLAGS_MINSIZEREL
		CMAKE_CXX_FLAGS_RELWITHDEBINFO)
	   if (${flag_var} MATCHES "/MD")
		  string (REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
	   endif ()
	endforeach ()

	# finds googletest built for the static runtime
	set (GTEST_MSVC_SEARCH "MT")
endif ()

# enables C++11 on GCC like compilers
if (CMAKE_COMPILER_IS_GNUCXX OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif ()

# finds Google Test
find_package (GTest)

# turn on if you want to generate tests
# turned on by default as long as GTest is found
option (GENERATE_TESTS "Generate Google Test tests" ${GTEST_FOUND})

# finds doxygen
find_package (Doxygen)

# turn on if you want to generate documentation
# turned on by default as long as Doxygen is found
option (GENERATE_DOCUMENTATION
	"Generate Doxygen documentation" ${DOXYGEN_FOUND})

# configures the project header
configure_file (src/singular/singular.h.in src/singular/singular.h @ONLY)

include_directories (
	src
	${PROJECT_BINARY_DIR}/src
	${GTEST_INCLUDE_DIRS})

add_executable (svd-test test/test.cpp)

# generates tests if necessary
if (GENERATE_TESTS)
	enable_testing ()

	add_executable (singular-test
		test/VectorTest.cpp
		test/MatrixTest.cpp
		test/DiagonalMatrixTest.cpp
		test/ReflectorTest.cpp
		test/RotatorTest.cpp
		test/SvdTest.cpp)

	# old Visual Studio needs a tweak
	if (MSVC AND MSVC_VERSION LESS 1800)
		set_target_properties (singular-test
			PROPERTIES COMPILE_FLAGS "-D_VARIADIC_MAX=10")
	endif ()

	target_link_libraries (singular-test ${GTEST_BOTH_LIBRARIES})

	add_test (singular-test-all singular-test)
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

# installs headers
install (FILES
	src/singular/DiagonalMatrix.h
	src/singular/Matrix.h
	src/singular/Reflector.h
	src/singular/Rotator.h
	src/singular/Svd.h
	src/singular/Vector.h
	${PROJECT_BINARY_DIR}/src/singular/singular.h
	DESTINATION include/singular)
