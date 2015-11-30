#include "Eigen/SVD"
#include "singular/Svd.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

/** Number of rows in an input matrix. */
static const int M = 60;

/** Number of columns in an input matrix. */
static const int N = 15;

/** Default number of iterations. */
static const int DEFAULT_ITERATION_COUNT = 1000;

/** Allowed error ratio. */
static const double ROUNDED_ERROR = 1.0e-8;

/** Matrix type for Eigen. */
typedef Eigen::Matrix< double, M, N > EigenMatrix;

/** SVD type for Eigen. */
typedef Eigen::JacobiSVD< EigenMatrix > EigenSVD;

/**
 * Runs Eigen's SVD over a matrix filled with given elements.
 *
 * @param elements
 *     Elements of the M x N matrix to be decomposed.
 */
static void runEigenSvd(const double elements[]) {
	// initializes the matrix
	EigenMatrix m;
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			m(i, j) = elements[i * N + j];
		}
	}
	// does decomposition
	EigenSVD svd(m, Eigen::ComputeFullU | Eigen::ComputeFullV);
}

/**
 * Runs singular's SVD over a matrix filled with given elements.
 *
 * @param elements
 *     Elements of the M x N matrix to be decomposed.
 */
static void runSingularSvd(const double elements[]) {
	// initializes the matrix
	singular::Matrix< M, N > m;
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			m(i, j) = elements[i * N + j];
		}
	}
	// does decomposition
	singular::Svd< M, N > svd;
	singular::Svd< M, N >::USV usv = svd.decomposeUSV(m);
}

/**
 * Runs a benchmark on a given decomposition function.
 *
 * @param numIterations
 *     Number of the iterations.
 * @param seed
 *     Seed for the random number generator.
 * @param f
 *     Decomposition function which takes an array `e` of elements to be
 *     decomposed, where the element at the ith row and jth column is given by
 *     `e[i * N + j]`
 */
double runBenchmark(
	int numIterations, unsigned int seed, void (*f)(const double[]))
{
	std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
	double elements[M * N];
	std::default_random_engine rnd(seed);
	std::uniform_real_distribution< double > dist(0.0, 1.0);
	for (int n = 0; n < numIterations; ++n) {
		// generates random elements
		std::generate_n(elements, M * N, [&rnd, &dist]() {
			return dist(rnd);
		});
		// runs the function
		f(elements);
	}
	std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
	std::chrono::duration< double > elapsed =
		std::chrono::duration_cast< std::chrono::duration< double > >(t2 - t1);
	return elapsed.count();
}

/**
 * Verifies results.
 *
 * @param numIterations
 *     Number of the iterations.
 * @param seed
 *     Seed for the random number generator.
 * @return
 *     Whether results of Eigen and singular match.
 */
bool verifyResults(int numIterations, unsigned int seed) {
	int numEigenDecompositionErrors = 0;
	double eigenDecompositionErrorSum = 0;
	int numEigenOrthonormalUErrors = 0;
	double eigenOrthonormalUErrorSum = 0;
	int numEigenOrthonormalVErrors = 0;
	double eigenOrthonormalVErrorSum = 0;
	int numDecompositionErrors = 0;
	double decompositionErrorSum = 0;
	int numOrthonormalUErrors = 0;
	double orthonormalUErrorSum = 0;
	int numOrthonormalVErrors = 0;
	double orthonormalVErrorSum = 0;
	int numSDiscrepancies = 0;
	double sErrorSum = 0;
	int numUDiscrepancies = 0;
	double uErrorSum = 0;
	int numVDiscrepancies = 0;
	double vErrorSum = 0;
	std::default_random_engine rnd(seed);
	std::uniform_real_distribution< double > dist(0.0, 1.0);
	double elements[M * N];
	for (int n = 0; n < numIterations; ++n) {
		// generates random elements
		std::generate_n(elements, M * N, [&rnd, &dist]() {
			return dist(rnd);
		});
		// does decomposition by Eigen
		Eigen::Matrix< double, M, N > m1;
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < N; ++j) {
				m1(i, j) = elements[i * N + j];
			}
		}
		EigenSVD svd1(m1, Eigen::ComputeFullU | Eigen::ComputeFullV);
		const EigenSVD::MatrixUType& u1 = svd1.matrixU();
		const EigenSVD::SingularValuesType& s1 = svd1.singularValues();
		const EigenSVD::MatrixVType& v1 = svd1.matrixV();
		// tests whether A = USV^T is satisfied
		{
			EigenMatrix sm;
			for (int i = 0; i < std::min(M, N); ++i) {
				sm(i, i) = s1(i);
			}
			EigenMatrix m3 = u1 * sm * v1.transpose();
			double mx = *std::max_element(elements, elements + M * N);
			for (int i = 0; i < M; ++i) {
				for (int j = 0; j < N; ++j) {
					double e = std::abs(m1(i, j) - m3(i, j));
					eigenDecompositionErrorSum += e;
					if (e / mx >= ROUNDED_ERROR) {
						++numEigenDecompositionErrors;
					}
				}
			}
		}
		// tests whether U*U^T = I is satisfied
		{
			Eigen::Matrix< double, M, M > eye = u1 * u1.transpose();
			for (int i = 0; i < M; ++i) {
				for (int j = 0; j < M; ++j) {
					double e;
					if (i == j) {
						e = std::abs(1.0 - eye(i, j));
					} else {
						e = std::abs(eye(i, j));
					}
					eigenOrthonormalUErrorSum += e;
					if (e >= ROUNDED_ERROR) {
						++numEigenOrthonormalUErrors;
					}
				}
			}
		}
		// tests whether V*V^T = I is satisfied
		{
			Eigen::Matrix< double, N, N > eye = v1 * v1.transpose();
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) {
					double e;
					if (i == j) {
						e = std::abs(1.0 - eye(i, j));
					} else {
						e = std::abs(eye(i, j));
					}
					eigenOrthonormalVErrorSum += e;
					if (e >= ROUNDED_ERROR) {
						++numEigenOrthonormalVErrors;
					}
				}
			}
		}
		// does decomposition by singular
		singular::Matrix< M, N > m2;
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < N; ++j) {
				m2(i, j) = elements[i * N + j];
			}
		}
		singular::Svd< M, N > svd2;
		singular::Svd< M, N >::USV usv = svd2.decomposeUSV(m2);
		const singular::Matrix< M, M >& u2 = singular::Svd< M, N >::getU(usv);
		const singular::DiagonalMatrix< M, N >& s2 =
			singular::Svd< M, N >::getS(usv);
		const singular::Matrix< N, N >& v2 = singular::Svd< M, N >::getV(usv);
		// tests whether A = USV^T is satisfied
		{
			double mx = *std::max_element(elements, elements + M * N);
			singular::Matrix< M, N > m3 = u2 * s2 * v2.transpose();
			for (int i = 0; i < M; ++i) {
				for (int j = 0; j < N; ++j) {
					double e = std::abs(m2(i, j) - m3(i, j));
					decompositionErrorSum += e;
					if (e / mx >= ROUNDED_ERROR) {
						++numDecompositionErrors;
					}
				}
			}
		}
		// tests whether U*U^T = I is satisfied
		{
			singular::Matrix< M, M > eye = u2 * u2.transpose();
			for (int i = 0; i < M; ++i) {
				for (int j = 0; j < M; ++j) {
					double e;
					if (i == j) {
						e = std::abs(1.0 - eye(i, j));
					} else {
						e = std::abs(eye(i, j));
					}
					orthonormalUErrorSum += e;
					if (e >= ROUNDED_ERROR) {
						++numOrthonormalUErrors;
					}
				}
			}
		}
		// tests whether V*V^T = I is satisfied
		{
			singular::Matrix< N, N > eye = v2 * v2.transpose();
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) {
					double e;
					if (i == j) {
						e = std::abs(1.0 - eye(i, j));
					} else {
						e = std::abs(eye(i, j));
					}
					orthonormalVErrorSum += e;
					if (e >= ROUNDED_ERROR) {
						++numOrthonormalVErrors;
					}
				}
			}
		}
		// compares singular values
		for (int i = 0; i < std::min(M, N); ++i) {
			double mx = std::max(std::abs(s1(i)), std::abs(s2(i, i)));
			double e = std::abs(s1(i) - s2(i, i));
			sErrorSum += e;
			if (mx != 0) {
				if (e / mx >= ROUNDED_ERROR) {
					++numSDiscrepancies;
				}
			}
		}
		// ignores sign flips
		std::vector< double > signs;
		for (int i = 0; i < std::min(M, N); ++i) {
			double sign = 1.0;
			if (u1(0, i) * u2(0, i) < 0.0) {
				sign = -1.0;
			}
			signs.push_back(sign);
		}
		// compares U matrices
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < std::min(M, N); ++j) {
				double mx = std::max(std::abs(u1(i, j)), std::abs(u2(i, j)));
				double e = std::abs(u1(i, j) - signs[j] * u2(i, j));
				uErrorSum += e;
				if (mx != 0) {
					if (e / mx >= ROUNDED_ERROR) {
						++numUDiscrepancies;
					}
				}
			}
		}
		// compares V matrices
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < std::min(M, N); ++j) {
				double mx = std::max(std::abs(v1(i, j)), std::abs(v2(i, j)));
				double e = std::abs(v1(i, j) - signs[j] * v2(i, j));
				vErrorSum += e;
				if (mx != 0) {
					if (e / mx >= ROUNDED_ERROR) {
						++numVDiscrepancies;
					}
				}
			}
		}
	}
	// shows statistics
	int numTotalElements = M * N * numIterations;
	int numTotalUElements = M * M * numIterations;
	int numTotalSElements = std::min(M, N) * numIterations;
	int numTotalVElements = N * N * numIterations;
	std::cout << "# of A=USV^T errors (Eigen): "
		<< numEigenDecompositionErrors << std::endl;
	std::cout << "mean A=USV^T error (Eigen): "
		<< (eigenDecompositionErrorSum / numTotalElements) << std::endl;
	std::cout << "# of A=USV^T errors (singular): "
		<< numDecompositionErrors << std::endl;
	std::cout << "mean A=USV^T error (singular): "
		<< (decompositionErrorSum / numTotalElements) << std::endl;
	std::cout << "# of U*U^T=I errors (Eigen): "
		<< numEigenOrthonormalUErrors << std::endl;
	std::cout << "mean U*U^T=I error (Eigen): "
		<< (eigenOrthonormalUErrorSum / numTotalUElements) << std::endl;
	std::cout << "# of U*U^T=I errors: " << numOrthonormalUErrors << std::endl;
	std::cout << "mean U*U^T=I error: "
		<< (orthonormalUErrorSum / numTotalUElements) << std::endl;
	std::cout << "# of V*V^T=I errors (Eigen): "
		<< numEigenOrthonormalVErrors << std::endl;
	std::cout << "mean V*V^T=I error (Eigen): "
		<< (eigenOrthonormalVErrorSum / numTotalVElements) << std::endl;
	std::cout << "# of V*V^T=I errors: " << numOrthonormalVErrors << std::endl;
	std::cout << "mean V*V^T=I error: "
		<< (orthonormalVErrorSum / numTotalVElements) << std::endl;
	std::cout << "# of S discrepancies: " << numSDiscrepancies << std::endl;
	std::cout << "mean S error: "
		<< (sErrorSum / numTotalSElements) << std::endl;
	std::cout << "# of U discrepancies: " << numUDiscrepancies << std::endl;
	std::cout << "mean U error: "
		<< (uErrorSum / (M * std::min(M, N) * numIterations)) << std::endl;
	std::cout << "# of V discrepancies: " << numVDiscrepancies << std::endl;
	std::cout << "mean V error: "
		<< (vErrorSum / (N * std::min(M, N) * numIterations)) << std::endl;
	return numDecompositionErrors == 0
		&& numOrthonormalUErrors == 0
		&& numOrthonormalVErrors == 0
		&& numSDiscrepancies == 0
		&& numUDiscrepancies == 0
		&& numVDiscrepancies == 0;
}

/**
 * Runs a benchmark.
 *
 * @param argv
 *     `argv[1]` is an optional number of iterations.
 */
int main(int argc, char** argv) {
	int numIterations = DEFAULT_ITERATION_COUNT;
	if (argc >= 2) {
		numIterations = atoi(argv[1]);
		if (numIterations <= 0) {
			std::cerr << "number of iterations must be a positive integer"
				<< " but " << argv[1];
			return 1;
		}
	}
	std::cout << "# of iterations: " << numIterations << std::endl;
	std::cout << "rounded error: " << ROUNDED_ERROR << std::endl;
	unsigned int seed =
		std::chrono::system_clock::now().time_since_epoch().count();
	// verifies the results
	bool verified = verifyResults(numIterations, seed);
	// measures times
	double eigenT = 0.0;
	double singularT = 0.0;
	double t = runBenchmark(numIterations, seed, &runSingularSvd);
	singularT += t;
	std::cout << "singular[1]: " << t << " seconds" << std::endl;
	t = runBenchmark(numIterations, seed, &runEigenSvd);
	eigenT += t;
	std::cout << "Eigen[1]: " << t << " seconds" << std::endl;
	t = runBenchmark(numIterations, seed, &runEigenSvd);
	eigenT += t;
	std::cout << "Eigen[2]: " << t << " seconds" << std::endl;
	t = runBenchmark(numIterations, seed, &runSingularSvd);
	singularT += t;
	std::cout << "singular[2]: " << t << " seconds" << std::endl;
	t = runBenchmark(numIterations, seed, &runSingularSvd);
	singularT += t;
	std::cout << "singular[3]: " << t << " seconds" << std::endl;
	t = runBenchmark(numIterations, seed, &runEigenSvd);
	eigenT += t;
	std::cout << "Eigen[3]: " << t << " seconds" << std::endl;
	t = runBenchmark(numIterations, seed, &runEigenSvd);
	eigenT += t;
	std::cout << "Eigen[4]: " << t << " seconds" << std::endl;
	t = runBenchmark(numIterations, seed, &runSingularSvd);
	singularT += t;
	std::cout << "singular[4]: " << t << " seconds" << std::endl;
	// shows statistics
	std::cout << "Eigen average: " << (eigenT / 4.0) << " seconds" << std::endl;
	std::cout << "singular average: "
		<< (singularT / 4.0) << " seconds" << std::endl;
	return verified ? 0 : 1;
}
