#pragma once
#include <string>

namespace co {

	typedef struct my_complex {
		float real;
		float imag;

		my_complex() {
			real = 0;
			imag = 0;
		}

		my_complex(float r, float i) {
			real = r;
			imag = i;
		}

	} cd;

	my_complex add(cd ,cd);
	my_complex sub(cd, cd);
	my_complex mul(cd, float);
	my_complex mul(cd, cd);
	my_complex e(cd);
	std::string toString(cd);
}
