#include "singular/singular.h"
#include "singular/Svd.h"

#include <iostream>

/** Runs a test. */
int main(int argc, char** argv) {
	const int M = 5;
	const int N = 4;
	const double DATA[] = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		4, 8, 3, 5,
		6, 7, 2, 1,
		9, 1, 3, 6
	};
	singular::Matrix< M, N > a;
	a.fill(DATA);
	singular::Svd< M, N >::USV usv = singular::Svd< M, N >::decomposeUSV(a);
	std::cout << "A = USV*" << std::endl;
	std::cout << std::endl;
	std::cout << "A = " << a << std::endl;
	std::cout << std::endl;
	std::cout << "U = " << singular::Svd< M, N >::getU(usv) << std::endl;
	std::cout << std::endl;
	std::cout << "S = " << singular::Svd< M, N >::getS(usv) << std::endl;
	std::cout << std::endl;
	std::cout << "V = " << singular::Svd< M, N >::getV(usv) << std::endl;
	std::cout << std::endl;
	std::cout << "singular version: " << SINGULAR_VERSION << std::endl;
	return 0;
}
