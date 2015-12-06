#include "singular/Svd.h"

#ifdef ENABLE_ARMADILLO
#include "armadillo"
#endif

#ifdef ENABLE_EIGEN
#include "Eigen/SVD"
#endif

#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

/** Number of rows in an input matrix. */
static const int M = 60;

/** Number of columns in an input matrix. */
static const int N = 50;

/** Minimum element value. */
static const double MIN_VALUE = -10;

/** Maximum element value. */
static const double MAX_VALUE = 10;

/** Default number of iterations. */
static const int DEFAULT_ITERATION_COUNT = 100;

/** Allowed error ratio. */
static const double ROUNDED_ERROR = 1.0e-12;

/** SVD configuration for singular. */
struct SingularSvd {
	/** Input matrix type. */
	typedef singular::Matrix< M, N > Matrix;

	/** Type for left-singular vectors. */
	typedef singular::Matrix< M, M > UMatrix;

	/** Type for singular values. */
	typedef singular::DiagonalMatrix< M, N > SMatrix;

	/** Type for right-singular vectors. */
	typedef singular::Matrix< N, N > VMatrix;

	/** Results of SVD. */
	singular::Svd< M, N >::USV usv;

	/**
	 * Performs singular value decomposition over given elements.
	 *
	 * @param elements
	 *     Elements of the M x N matrix to be decomposed.
	 */
	void operator ()(const double elements[]) {
		// initializes the matrix
		singular::Matrix< M, N > m;
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < N; ++j) {
				m(i, j) = elements[i * N + j];
			}
		}
		// does decomposition
		singular::Svd< M, N > svd;
		this->usv = svd.decomposeUSV(m);
	}

	/** Returns the last computed left-singular matrix. */
	inline const UMatrix& getU() const {
		return singular::Svd< M, N >::getU(this->usv);
	}

	/** Returns the transposed last computed left-singular matrix. */
	inline UMatrix getUT() const {
		return singular::Svd< M, N >::getU(this->usv).transpose();
	}

	/** Returns the last computed singular values. */
	inline SMatrix getS() const {
		return singular::Svd< M, N >::getS(this->usv).clone();
	}

	/** Returns the last computed right-singular vectors. */
	inline const VMatrix& getV() const {
		return singular::Svd< M, N >::getV(this->usv);
	}

	/** Returns the transposed last computed right-singular vectors. */
	inline VMatrix getVT() const {
		return singular::Svd< M, N >::getV(this->usv).transpose();
	}
};

#ifdef ENABLE_EIGEN
/** SVD configuration for Eigen. */
struct EigenSvd {
	/** Input matrix type. */
	typedef Eigen::Matrix< double, M, N > Matrix;

	/** Type for left-singular vectors. */
	typedef Eigen::JacobiSVD< Matrix >::MatrixUType UMatrix;

	/** Type for singular values. */
	typedef Matrix SMatrix;

	/** Type for right-singular vectors. */
	typedef Eigen::JacobiSVD< Matrix >::MatrixVType VMatrix;

	/** Results of SVD. */
	Eigen::JacobiSVD< Matrix > svd;

	/**
	 * Performs singular value decomposition over given elements.
	 *
	 * @param elements
	 *     Elements of the M x N matrix to be decomposed.
	 */
	void operator ()(const double elements[]) {
		// initializes the matrix
		Matrix m;
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < N; ++j) {
				m(i, j) = elements[i * N + j];
			}
		}
		// does decomposition
		this->svd.compute(m, Eigen::ComputeFullU | Eigen::ComputeFullV);
	}

	/** Returns the last computed left-singular vectors. */
	inline const UMatrix& getU() const {
		return this->svd.matrixU();
	}

	/** Returns the transposed last computed left-singular vectors. */
	inline UMatrix getUT() const {
		return this->svd.matrixU().transpose();
	}

	/**
	 * Returns the last computed singular values.
	 * Converts singular value vector into an M x N diagonal matrix.
	 */
	SMatrix getS() const {
		const Eigen::JacobiSVD< Matrix >::SingularValuesType& ss =
			this->svd.singularValues();
		SMatrix sm;
		sm.setZero();
		for (int i = 0; i < std::min(M, N); ++i) {
			sm(i, i) = ss(i);
		}
		return sm;
	}

	/** Returns the last computed right-singular vectors. */
	inline const VMatrix& getV() const {
		return this->svd.matrixV();
	}

	/** Returns the transposed last computed right-singular vectors. */
	inline VMatrix getVT() const {
		return this->svd.matrixV().transpose();
	}
};
#endif

#ifdef ENABLE_ARMADILLO
/** SVD configuration for Armadillo. */
struct ArmadilloSvd {
	/** Input matrix type. */
	typedef arma::mat Matrix;

	/** Type for left-singular vectors. */
	typedef arma::mat UMatrix;

	/** Type for singular values. */
	typedef arma::mat SMatrix;

	/** Type for right-singular vectors. */
	typedef arma::mat VMatrix;

	/** Stores left-singular vectors. */
	UMatrix u;

	/** Stores singular values. */
	arma::vec s;

	/** Stores right-singular vectors. */
	VMatrix v;

	/**
	 * Performs singular value decomposition over given elements.
	 *
	 * @param elements
	 *     Elements of the M x N matrix to be decomposed.
	 */
	void operator ()(const double elements[]) {
		// initializes the matrix
		arma::mat m(M, N);
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < N; ++j) {
				m(i, j) = elements[i * N + j];
			}
		}
		// does decomposition
		arma::svd(this->u, this->s, this->v, m);
	}

	/** Returns the last computed left-singular vectors. */
	inline const UMatrix& getU() const {
		return this->u;
	}

	/** Returns the transposed last computed left-singular vectors. */
	inline UMatrix getUT() const {
		return this->u.t();
	}

	/**
	 * Returns the last computed singular values.
	 * Converts singular value vector into M x N diagonal matrix.
	 */
	SMatrix getS() const {
		SMatrix sm(M, N, arma::fill::zeros);
		for (int i = 0; i < std::min(M, N); ++i) {
			sm(i, i) = this->s(i);
		}
		return sm;
	}

	/** Returns the last computed right-singular vectors. */
	inline const VMatrix& getV() const {
		return this->v;
	}

	/** Returns the transposed last computed right-singular vectors. */
	inline VMatrix getVT() const {
		return this->v.t();
	}
};
#endif

/**
 * A benchmark function for a given algorithm.
 *
 * @tparam Algorithm
 *     Algorithm to be evaluated.
 */
template < typename Algorithm >
struct Benchmark {
	/** Number of iterations. */
	int numIterations;

	/** Seed for the random number generator. */
	int seed;

	/**
	 * Configures a benchmark.
	 *
	 * @param numIterations
	 *     Number of iterations.
	 * @param seed
	 *     Seed for the random number generator.
	 */
	Benchmark(int numIterations, int seed)
		: numIterations(numIterations), seed(seed) {}

	/** Runs the given algorithm several times. */
	void operator ()() const {
		Algorithm algo;
		double elements[M * N];
		std::default_random_engine rnd(this->seed);
		std::uniform_real_distribution< double > dist(MIN_VALUE, MAX_VALUE);
		for (int n = 0; n < this->numIterations; ++n) {
			// generates random elements
			std::generate_n(elements, M * N, [&rnd, &dist]() {
				return dist(rnd);
			});
			// runs the algorithm
			algo(elements);
		}
	}
};

/** Stopwatch to evaluate an algorithm. */
class Stopwatch {
private:
	/** Clock type used to measure elapsed time. */
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ == 4 && __GNUC_MINOR__ < 7
	// GCC prior to 4.7 does not have steady_clock
	typedef std::chrono::monotonic_clock ClockType;
#else
	typedef std::chrono::steady_clock ClockType;
#endif

	/** Time point when this stopwatch has started. */
	ClockType::time_point sT;

	/** Measured lap times. */
	std::vector< double > lapTimes;
public:
	/** Initializes a stopwatch. */
	Stopwatch() {
		this->lapTimes.reserve(10);
	}

	/**
	 * Measures the time to run a given function.
	 *
	 * @tparam Fn
	 *     Type of the function to run.
	 * @param fn
	 *     Function or function like object to run.
	 */
	template < typename Fn >
	void measure(Fn& fn) {
		this->start();
		fn();
		this->stop();
	}

	/** Prints the statistics on the standard output. */
	void printStatistics() const {
		double sum = 0.0;
		for (size_t i = 0; i < this->lapTimes.size(); ++i) {
			std::cout << "lap time[" << i << "]: "
				<< this->lapTimes[i] << " seconds" << std::endl;
			sum += this->lapTimes[i];
		}
		std::cout << "mean lap time: "
			<< (sum / this->lapTimes.size()) << " seconds" << std::endl;
	}
private:
	/** Starts measuring time. */
	void start() {
		this->sT = std::chrono::steady_clock::now();
	}

	/**
	 * Stops measuring time and records the elapsed time since this stopwatch
	 * has last started.
	 */
	void stop() {
		ClockType::time_point eT = ClockType::now();
		std::chrono::duration< double > elapsed =
			std::chrono::duration_cast< std::chrono::duration< double > >(eT - this->sT);
		this->lapTimes.push_back(elapsed.count());
	}
};

/**
 * Verifier for a given SVD algorithm.
 *
 * @tparam Algorithm
 *     SVD algorithm to be verified.
 *     Must have a default constructor.
 */
template < typename Algorithm >
class SvdVerifier {
private:
	/** Algorithm to be verified. */
	Algorithm algo;

	/** Number of verifications. */
	int numVerifications;

	/** Number of reconstruction errors. */
	int numReconstructionErrors;

	/** Sum of reconstruction errors. */
	double reconstructionErrorSum;

	/** Number of orthonormal test errors over left-singular vectors. */
	int numOrthonormalUErrors;

	/** Sum of orthonormal test errors over left-singular vectors. */
	double orthonormalUErrorSum;

	/** Number of orthonormal test errors over right-singular vectors. */
	int numOrthonormalVErrors;

	/** Sum of orthonormal test errors over right-singular vectors. */
	double orthonormalVErrorSum;

	/** Number of reference singular value discrepancies. */
	int numSingularValueDiscrepancies;

	/** Sum of singular value errors. */
	double singularValueErrorSum;
public:
	/** Initializes. */
	SvdVerifier()
		: numVerifications(0),
		  numReconstructionErrors(0),
		  reconstructionErrorSum(0),
		  numOrthonormalUErrors(0),
		  orthonormalUErrorSum(0),
		  numOrthonormalVErrors(0),
		  orthonormalVErrorSum(0),
		  numSingularValueDiscrepancies(0),
		  singularValueErrorSum(0) {}

	/** Returns the algorithm object being verified. */
	inline Algorithm& getAlgorithm() {
		return this->algo;
	}

	/**
	 * Verifies the algorithm with given elements.
	 *
	 * @param elements
	 *     Elements of the M x N matrix.
	 */
	void verify(const double elements[]) {
		++this->numVerifications;
		this->algo(elements);
		// tests reconstruction A = USV*
		{
			double mx = *std::max_element(
				elements, elements + M * N, [](double a, double b) {
					return std::abs(a) < std::abs(b);
				});
			if (mx == 0.0) {
				mx = 1.0;
			}
			typename Algorithm::Matrix a =
				this->algo.getU() * this->algo.getS() * this->algo.getVT();
			for (int i = 0; i < M; ++i) {
				for (int j = 0; j < N; ++j) {
					double ref = elements[i * N + j];
					double e = std::abs(a(i, j) - ref);
					this->reconstructionErrorSum += e;
					if (e / mx >= ROUNDED_ERROR) {
						++this->numReconstructionErrors;
					}
				}
			}
		}
		// tests UU* = I
		{
			typename Algorithm::UMatrix eye =
				this->algo.getU() * this->algo.getUT();
			for (int i = 0; i < M; ++i) {
				for (int j = 0; j < M; ++j) {
					double e;
					if (i == j) {
						e = std::abs(eye(i, j) - 1.0);
					} else {
						e = std::abs(eye(i, j));
					}
					this->orthonormalUErrorSum += e;
					if (e >= ROUNDED_ERROR) {
						++this->numOrthonormalUErrors;
					}
				}
			}
		}
		// tests VV* = I
		{
			typename Algorithm::VMatrix eye =
				this->algo.getV() * this->algo.getVT();
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) {
					double e;
					if (i == j) {
						e = std::abs(eye(i, j) - 1.0);
					} else {
						e = std::abs(eye(i, j));
					}
					this->orthonormalVErrorSum += e;
					if (e >= ROUNDED_ERROR) {
						++this->numOrthonormalVErrors;
					}
				}
			}
		}
	}

	/**
	 * Compares singular values with given reference values.
	 *
	 * This must follow the call of the `verify` function.
	 *
	 * @tparam Matrix
	 *     Matrix type of the reference singular values.
	 * @param ref
	 *     Reference singular values.
	 */
	template < typename Matrix >
	void compareSingularValues(const Matrix& ref) {
		typename Algorithm::SMatrix s = this->algo.getS();
		for (int i = 0; i < std::min(M, N); ++i) {
			double mx = std::max(std::abs(s(i, i)), std::abs(ref(i, i)));
			double e = std::abs(s(i, i) - ref(i, i));
			this->singularValueErrorSum += e;
			if (e / mx >= ROUNDED_ERROR) {
				++this->numSingularValueDiscrepancies;
			}
		}
	}

	/**
	 * Returns whether the algorithm has been verified.
	 *
	 * The algorithm is considered verified if it passes all of the following
	 * tests,
	 *  - No reconstruction errors
	 *  - No orthonormal left-singular vector errors
	 *  - No orthonormal right-singular vector errors
	 *
	 * @return
	 *     Whether the algorithm has been verified.
	 */
	inline bool isVerified() const {
		return this->numReconstructionErrors == 0
			&& this->numOrthonormalUErrors == 0
			&& this->numOrthonormalVErrors == 0;
	}

	/** Prints statistics to the standard output. */
	void printStatistics() {
		double numReconstructionTests = M * N * this->numVerifications;
		double numOrthonormalUTests = M * M * this->numVerifications;
		double numOrthonormalVTests = N * N * this->numVerifications;
		double numSingularValueTests = std::min(M, N) * this->numVerifications;
		std::cout
			<< "# of reconstruction errors: "
			<< this->numReconstructionErrors
			<< "  mean error: "
			<< (this->reconstructionErrorSum / numReconstructionTests)
			<< std::endl;
		std::cout
			<< "# of orthonormal U errors: "
			<< this->numOrthonormalUErrors
			<< "  mean error: "
			<< (this->orthonormalUErrorSum / numOrthonormalUTests)
			<< std::endl;
		std::cout
			<< "# of orthonormal V errors: "
			<< this->numOrthonormalVErrors
			<< "  mean error: "
			<< (this->orthonormalVErrorSum / numOrthonormalVTests)
			<< std::endl;
		std::cout
			<< "# of singular value discrepancies: "
			<< this->numSingularValueDiscrepancies
			<< "  mean error: "
			<< (this->singularValueErrorSum / numSingularValueTests)
			<< std::endl;
	}
};

/**
 * Verifies results.
 *
 * @param numIterations
 *     Number of the iterations.
 * @param seed
 *     Seed for the random number generator.
 * @return
 *     Whether singular's results are verified.
 */
bool verifyResults(int numIterations, unsigned int seed) {
	SvdVerifier< SingularSvd > singularVerifier;
#ifdef ENABLE_EIGEN
	SvdVerifier< EigenSvd > eigenVerifier;
#endif
#ifdef ENABLE_ARMADILLO
	SvdVerifier< ArmadilloSvd > armadilloVerifier;
#endif
	std::default_random_engine rnd(seed);
	std::uniform_real_distribution< double > dist(MIN_VALUE, MAX_VALUE);
	double elements[M * N];
	for (int n = 0; n < numIterations; ++n) {
		std::generate_n(elements, M * N, [&rnd, &dist]() {
			return dist(rnd);
		});
		singularVerifier.verify(elements);
#ifdef ENABLE_EIGEN
		eigenVerifier.verify(elements);
		eigenVerifier.compareSingularValues(
			singularVerifier.getAlgorithm().getS());
#endif
#ifdef ENABLE_ARMADILLO
		armadilloVerifier.verify(elements);
		armadilloVerifier.compareSingularValues(
			singularVerifier.getAlgorithm().getS());
#endif
	}
	// shows statistics
	std::cout << "singular" << std::endl;
	singularVerifier.printStatistics();
	std::cout << std::endl;
#ifdef ENABLE_EIGEN
	std::cout << "Eigen" << std::endl;
	eigenVerifier.printStatistics();
	std::cout << std::endl;
#endif
#ifdef ENABLE_ARMADILLO
	std::cout << "Armadillo" << std::endl;
	armadilloVerifier.printStatistics();
	std::cout << std::endl;
#endif
	return singularVerifier.isVerified();
}

/**
 * Runs a benchmark.
 *
 * @param argv
 *     `argv[1]` is an optional number of iterations.
 */
int main(int argc, char** argv) {
	// gets the number of iterations if given
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
	std::cout << "min value: " << MIN_VALUE << std::endl;
	std::cout << "max value: " << MAX_VALUE << std::endl;
	std::cout << std::endl;
	std::cout << "verifying results ..." << std::endl;
	std::cout << std::endl;
	unsigned int seed =
		std::chrono::system_clock::now().time_since_epoch().count();
	// verifies the results
	bool verified = verifyResults(numIterations, seed);
	// runs benchmarks
	Benchmark< SingularSvd > singularBenchmark(numIterations, seed);
	Stopwatch singularWatch;
#ifdef ENABLE_EIGEN
	Benchmark< EigenSvd > eigenBenchmark(numIterations, seed);
	Stopwatch eigenWatch;
#endif
#ifdef ENABLE_ARMADILLO
	Benchmark< ArmadilloSvd > armadilloBenchmark(numIterations, seed);
	Stopwatch armadilloWatch;
#endif
	std::cout << "measuring processing time ..." << std::endl;
	// round 1
	std::cout << "round 1" << std::endl;
	singularWatch.measure(singularBenchmark);
#ifdef ENABLE_EIGEN
	eigenWatch.measure(eigenBenchmark);
#endif
#ifdef ENABLE_ARMADILLO
	armadilloWatch.measure(armadilloBenchmark);
#endif
	// round 2
	std::cout << "round 2" << std::endl;
	singularWatch.measure(singularBenchmark);
#ifdef ENABLE_ARMADILLO
	armadilloWatch.measure(armadilloBenchmark);
#endif
#ifdef ENABLE_EIGEN
	eigenWatch.measure(eigenBenchmark);
#endif
	// round 3
	std::cout << "round 3" << std::endl;
#ifdef ENABLE_EIGEN
	eigenWatch.measure(eigenBenchmark);
#endif
	singularWatch.measure(singularBenchmark);
#ifdef ENABLE_ARMADILLO
	armadilloWatch.measure(armadilloBenchmark);
#endif
	// round 4
	std::cout << "round 4" << std::endl;
#ifdef ENABLE_EIGEN
	eigenWatch.measure(eigenBenchmark);
#endif
#ifdef ENABLE_ARMADILLO
	armadilloWatch.measure(armadilloBenchmark);
#endif
	singularWatch.measure(singularBenchmark);
	// round 5
	std::cout << "round 5" << std::endl;
#ifdef ENABLE_ARMADILLO
	armadilloWatch.measure(armadilloBenchmark);
#endif
#ifdef ENABLE_EIGEN
	eigenWatch.measure(eigenBenchmark);
#endif
	singularWatch.measure(singularBenchmark);
	// round 6
	std::cout << "round 6" << std::endl;
#ifdef ENABLE_ARMADILLO
	armadilloWatch.measure(armadilloBenchmark);
#endif
	singularWatch.measure(singularBenchmark);
#ifdef ENABLE_EIGEN
	eigenWatch.measure(eigenBenchmark);
#endif
	std::cout << std::endl;
	// shows statistics
	std::cout << "singular: " << std::endl;
	singularWatch.printStatistics();
	std::cout << std::endl;
#ifdef ENABLE_EIGEN
	std::cout << "Eigen: " << std::endl;
	eigenWatch.printStatistics();
	std::cout << std::endl;
#endif
#ifdef ENABLE_ARMADILLO
	std::cout << "Armadillo: " << std::endl;
	armadilloWatch.printStatistics();
	std::cout << std::endl;
#endif
	return verified ? 0 : 1;
}
