cmake_minimum_required (VERSION 2.8.4)

project (singular)

set (PROJECT_MAJOR 0)
set (PROJECT_MINOR 1)
set (PROJECT_PATCH 3-SNAPSHOT)
set (PROJECT_VERSION "${PROJECT_MAJOR}.${PROJECT_MINOR}.${PROJECT_PATCH}")

# turn off if you do not want to compile benchmark
option (ENABLE_BENCHMARK "Enable benchmark" on)

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
if (CMAKE_COMPILER_IS_GNUCXX)
	# checks the GCC version
	#  - 4.7 or higher needs -std=c++11
	#  - 4.6 or lower needs -std=c++0x
	# http://stackoverflow.com/questions/4058565/check-gcc-minor-in-cmake
	execute_process (
		COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
		OUTPUT_VARIABLE GCC_VERSION)
	if (GCC_VERSION VERSION_EQUAL "4.7" OR GCC_VERSION VERSION_GREATER "4.7")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	else ()
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
	endif ()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
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

# benchmark
if (ENABLE_BENCHMARK)
	add_executable (benchmark test/benchmark.cpp)

	# additional compiler flags for benchmark
	set (BENCHMARK_COMPILE_FLAGS "")

	# enables benchmarks on Eigen if EIGEN_INCLUDE is set
	if (EIGEN_INCLUDE)
		if (CMAKE_VERSION VERSION_EQUAL "3.0"
				OR CMAKE_VERSION VERSION_GREATER "3.0")
			target_include_directories (benchmark
				PUBLIC ${EIGEN_INCLUDE})
		else ()
			include_directories (${EIGEN_INCLUDE})
		endif ()
		set (BENCHMARK_COMPILE_FLAGS
			"${BENCHMARK_COMPILE_FLAGS} -DENABLE_EIGEN")
	endif ()

	# enables benchmarks on Armadillo if ARMADILLO_ROOT is set
	if (ARMADILLO_ROOT)
		# locates the Armadillo, LAPACK and BLAS libraries
		find_library (ARMADILLO_LIB armadillo HINTS ${ARMADILLO_ROOT}/lib)
		if (APPLE)
			# uses Accelerate framework on Apple products
			set (LAPACK_FOUND true)
			set (BLAS_FOUND true)
			set (LAPACK_LIBRARIES "-framework Accelerate")
			set (BLAS_LIBRARIES "")
		else ()
			find_package (LAPACK)
			find_package (BLAS)
		endif ()
		if (ARMADILLO_LIB AND LAPACK_FOUND AND BLAS_FOUND)
			if (CMAKE_VERSION VERSION_EQUAL "3.0"
					OR CMAKE_VERSION VERSION_GREATER "3.0")
				target_include_directories (benchmark
					PUBLIC ${ARMADILLO_ROOT}/include)
			else ()
				include_directories (${ARMADILLO_ROOT}/include)
			endif ()
			target_link_libraries (benchmark
				${ARMADILLO_LIB}
				${LAPACK_LIBRARIES}
				${BLAS_LIBRARIES})
			set (BENCHMARK_COMPILE_FLAGS
				"${BENCHMARK_COMPILE_FLAGS} -DENABLE_ARMADILLO")
		else ()
			if (NOT ARMADILLO_LIB_FOUND)
				message (WARNING
					"no Armadillo library found in ${ARMADILLO_ROOT}/lib")
			endif ()
			if (NOT LAPACK_FOUND)
				message (WARNING "no LAPACK found")
			endif ()
			if (NOT BLAS_FOUND)
				message (WARNING "no BLAS found")
			endif ()
		endif ()
	endif ()

	if (BENCHMARK_COMPILE_FLAGS)
		set_target_properties (benchmark
			PROPERTIES COMPILE_FLAGS ${BENCHMARK_COMPILE_FLAGS})
	endif ()
endif ()

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
