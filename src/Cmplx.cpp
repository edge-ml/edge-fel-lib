#include "Cmplx.h"
#include <cmath>

using namespace co;
using namespace std;

my_complex co::add(cd a, cd b) {
	cd res{};
	res.real = a.real + b.real;
	res.imag = a.imag + b.imag;
	return res;
}

my_complex co::sub(cd a, cd b) {
	cd res{};
	res.real = a.real - b.real;
	res.imag = a.imag - b.imag;
	return res;
}

my_complex co::mul(cd a, float b) {
	cd res{};
	res.real = a.real * b;
	res.imag = a.imag * b;
	return res;
}

my_complex co::mul(cd a, cd b) {
	cd res{};
	res.real = a.real * b.real - a.imag * b.imag;
	res.imag = a.real * b.imag + a.imag * b.real;
	return res;
}

my_complex co::e(cd a) {
	cd res{};
	res.real = exp(a.real) * cos(a.imag);
	res.imag = exp(a.real) * sin(a.imag);
	return res;
}

string co::toString(cd a) {
	return "(" + to_string(a.real) + ", " + to_string(a.imag) + ")";
}
