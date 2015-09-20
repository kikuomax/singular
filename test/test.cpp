#include "singular/svd.h"

#include <iostream>

/** Runs a test. */
int main(int argc, char** argv) {
	singular::Matrix< 5, 4 > a;
	const double DATA[] = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		4, 8, 3, 5,
		6, 7, 2, 1,
		9, 1, 3, 6
	};
	a.fill(DATA);
	std::tuple< singular::Matrix< 5, 5 >,
				singular::Matrix< 5, 4 >,
				singular::Matrix< 4, 4 > > usv = singular::svdUSV(a);
	std::cout << "A = USV*" << std::endl;
	std::cout << std::endl;
	std::cout << "A = " << a << std::endl;
	std::cout << std::endl;
	std::cout << "U = " << std::get< 0 >(usv) << std::endl;
	std::cout << std::endl;
	std::cout << "S = " << std::get< 1 >(usv) << std::endl;
	std::cout << std::endl;
	std::cout << "V = " << std::get< 2 >(usv) << std::endl;
	return 0;
}
