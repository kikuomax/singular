*singular* is a C++ library which provides functions and classes to perform [singular value decomposition (SVD)](https://en.wikipedia.org/wiki/Singular_value_decomposition).

[![Build Status](https://travis-ci.org/kikuomax/singular.svg?branch=master)](https://travis-ci.org/kikuomax/singular)

The goals of this reinvention are
 - Learning Francis's algorithm
 - Learning modern C++

Prerequisites
=============

You need the following software installed,
 - [CMake](https://cmake.org) 2.8.4 or higher
 - C++ compiler

Getting started
===============

*singular* consists of only header files.
The [CMake script](/CMakeLists.txt) does everything you need to install *singular*.

Please take the following steps,

 1. Clone the repository anywhere you like and move down to it.

	```shell
	git clone https://github.com/kikuomax/singular.git
	cd singular
	```

 2. Create a `build` directory and move down to it.

	```shell
	mkdir build
	cd build
	```

 3. Configure the project.
    Suppose you want to install *singular* into the directory `install-path`.

	```shell
	cmake -DCMAKE_INSTALL_PREFIX=install-path ..
	```

 4. Install headers.

    ```shell
	cmake --build . --target install
	```

 5. You will find the headers in the following directory,
     - `install-path/include/singular`

Creation of the `build` directory (step 2) is not necessary, but it prevents this directory being messy.

An [example program](/test/test.cpp) shows how to use *singular*.

Testing singular
================

If you have [googletest](https://code.google.com/p/googletest/) installed, unittests for *singular* are built by the CMake script as well.
You can run unittests by the following command after running the build step of *singular*,

```shell
ctest -V
```

pthread related errors
----------------------

With the following configuration, I got linker errors related to pthread.
 - Ubuntu 14.04.3
 - gcc 4.8.4
 - googletest 1.7.0

One solution is to add a `-pthread` compiler flag followed by other objects and libraries.
[This conversation](https://github.com/google/googletest/issues/391) may help you.
You can specify the `-pthread` flag at the configuration step like,

```shell
cmake -DCMAKE_CXX_FLAGS=-pthread ..
```

Another solution is to compile googletest without pthread.
This can be done by turning on the `gtest_disable_pthreads` option at the configuration of **googletest**.

```shell
cmake -Dgtest_disable_pthreads=on .
```

Benchmarking
============

[`test/benchmark.cpp`](/test/benchmark.cpp) is a simple benchmarking program.
This program measures time needed to perform SVD over randomly generated matrices.
It also evaluates the following qualities,
 - Reversibility (`A = USV*`)
 - Orthonormality of left-singular vectors (`UU* = I`)
 - Orthonormality of right-singular vectors (`VV* = I`)

If you do not want to compile the benchmark, turn off the `ENABLE_BENCHMARK` option.

[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [Armadillo](http://arma.sourceforge.net) can optionally be evaluated.
This feature is initially disabled.

To enable the benchmark on Eigen, please set the `EIGEN_INCLUDE` option to the directory where header files of Eigen are installed.

To enable the benchmark on Armadillo, please set the `ARMADILLO_ROOT` option to the directory where Armadillo is installed.
Armadillo needs [LAPACK](http://www.netlib.org/lapack/) and [BLAS](http://www.netlib.org/blas/) installed.
On Mac OS X, you can link LAPACK and BLAS by adding a compiler flag `-framework Accelerate`.

The following are the results on my PC (MacBook Pro, 2.4GHz Intel Core i5).

```
# of rows: 60
# of columns: 50
# of iterations: 100
rounded error: 1e-12
min value: -10
max value: 10
singular version: 0.1.3
Eigen version: 3.2.7
Armadillo version: 6.200.5 (Midnight Blue Deluxe)

verifying results ...

singular
# of reconstruction errors: 0  mean error: 1.96851e-14
# of orthonormal U errors: 0  mean error: 2.04803e-16
# of orthonormal V errors: 0  mean error: 2.42038e-16
# of singular value discrepancies: 0  mean error: 0

Eigen
# of reconstruction errors: 0  mean error: 1.13796e-13
# of orthonormal U errors: 0  mean error: 8.86169e-16
# of orthonormal V errors: 0  mean error: 4.77086e-16
# of singular value discrepancies: 0  mean error: 4.03472e-13

Armadillo
# of reconstruction errors: 0  mean error: 1.08298e-14
# of orthonormal U errors: 0  mean error: 1.60648e-16
# of orthonormal V errors: 0  mean error: 1.97627e-16
# of singular value discrepancies: 0  mean error: 5.28169e-14

measuring processing time ...
round 1
round 2
round 3
round 4
round 5
round 6

singular: 
lap time[0]: 3.09783 seconds
lap time[1]: 3.08843 seconds
lap time[2]: 3.08597 seconds
lap time[3]: 3.07436 seconds
lap time[4]: 3.07892 seconds
lap time[5]: 3.08201 seconds
mean lap time: 3.08459 seconds

Eigen: 
lap time[0]: 6.63677 seconds
lap time[1]: 6.61933 seconds
lap time[2]: 6.64333 seconds
lap time[3]: 6.66512 seconds
lap time[4]: 6.64167 seconds
lap time[5]: 6.62753 seconds
mean lap time: 6.63896 seconds

Armadillo: 
lap time[0]: 0.0938214 seconds
lap time[1]: 0.0937346 seconds
lap time[2]: 0.0936364 seconds
lap time[3]: 0.0936577 seconds
lap time[4]: 0.100055 seconds
lap time[5]: 0.10046 seconds
mean lap time: 0.0958942 seconds
```

Yes! Armadillo could be the best choice if you can use LAPACK and BLAS on your application.

Generating documentation
========================

If you have [Doxygen](http://www.doxygen.org) installed, documentation of *singular* is generated by the CMake script as well.
You will find the `api` directory created by CMake.

License
=======

[MIT License](https://opensource.org/licenses/MIT)

References
==========

Most of the code was derived from the following great book.

 - Watkins, David S. *Fundamentals of Matrix Computations*. 3rd ed. Hoboken: John Wiley & Sons, 2010


Test results of SVD were generated at the following site.

 - http://www.bluebit.gr/matrix-calculator/
