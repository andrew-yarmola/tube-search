#include "AJRC.h"
#include "types.hh"

template<>
void print_type<const AJRC>(const AJRC& x) {
	fprintf(stderr, "f: %.20f + %.20f I\n\
z0: %.20f + %.20f I   z1: %.20f + %.20f I\n\
z2: %.20f + %.20f I   z3: %.20f + %.20f I\n\
err: %.40f\n\
size: %.40f\n\
absLB: %.20f | hex %s\n\
abdUB: %.20f | hex %s\n", x.f.re, x.f.im,
		   x.z0.re, x.z0.im, x.z1.re, x.z1.im,
		   x.z2.re, x.z2.im, x.z3.re, x.z3.im, x.e, x.size,
       absLB(x), double_to_hex(absLB(x)).c_str(),
       absUB(x), double_to_hex(absUB(x)).c_str());
}

template<>
void print_type<AJRC>(AJRC& x) {
  print_type((const AJRC) x);
}

template<>
void print_center<const AJRC>(const AJRC& x) {
	printf("f: %f + %f I\nabsLB: %f, abdUB: %f\n", x.f.re, x.f.im, absLB(x), absUB(x));
}

template<>
bool sort_comp<AJRC>(const AJRC& a, const AJRC& b) {
  return a.e < b.e;
}

const AJRC operator*(const AJRC&x,const AJRC&y) {

	double xdist = size(x);
	double ydist = size(y);
	double ax = absUB(x.f), ay = absUB(y.f);
	AComplex r_f  = x.f * y.f;
	AComplex r_z0 = x.f * y.z0 + x.z0 * y.f;
	AComplex r_z1 = x.f * y.z1 + x.z1 * y.f;
	AComplex r_z2 = x.f * y.z2 + x.z2 * y.f;
	AComplex r_z3 = x.f * y.z3 + x.z3 * y.f;
	double A = (xdist + x.e) * (ydist + y.e);
	double B = ax * y.e + ay * x.e;
	double C = (r_f.e + (r_z0.e + r_z1.e)) + (r_z2.e + r_z3.e);
	double r_error = (1 + 3 * EPS) * ((A + B) + C);

    return AJRC(r_f.z, r_z0.z, r_z1.z, r_z2.z, r_z3.z, r_error);

}
const AJRC operator/(const AJRC&x,const AJRC&y) {

	double xdist = size(x);
	double ydist = size(y);
	double ax = absUB(x.f), ay = absLB(y.f);
	double D = ay - (1 + EPS) * (y.e + ydist);
	if(!(D > 0)) return AJRC(0, 0, 0, 0, 0, 0, 0, infinity());
	AComplex den = (y.f * y.f);
	AComplex r_f  = x.f / y.f;
	AComplex r_z0 = (x.z0 * y.f - x.f * y.z0) / den;
	AComplex r_z1 = (x.z1 * y.f - x.f * y.z1) / den;
	AComplex r_z2 = (x.z2 * y.f - x.f * y.z2) / den;
	AComplex r_z3 = (x.z3 * y.f - x.f * y.z3) / den;
	double A = (ax + (xdist + x.e)) / D;
	double B = (ax / ay + xdist / ay) + (ydist * ax) / (ay * ay);
	double C = (r_f.e + (r_z0.e + r_z1.e)) + (r_z2.e + r_z3.e);
	double r_error = (1 + 3 * EPS) * (((1 + 3 * EPS) * A - (1 - 3 * EPS) * B) + C);

    return AJRC(r_f.z, r_z0.z, r_z1.z, r_z2.z, r_z3.z, r_error);

}
const AJRC operator/(double x,const AJRC&y) {

	double ydist = size(y);
	double ax = fabs(x), ay = absLB(y.f);
	double D = ay - (1 + EPS) * (y.e + ydist);
	if(!(D > 0))return AJRC(0, 0, 0, 0, 0, 0, 0, infinity());
	AComplex den = (y.f * y.f);
	AComplex r_f =  x / y.f;
	AComplex r_z0 = (-x * y.z0) / den;
	AComplex r_z1 = (-x * y.z1) / den;
	AComplex r_z2 = (-x * y.z2) / den;
	AComplex r_z3 = (-x * y.z3) / den;
	double B = ax / ay + (ydist * ax) / ( ay * ay);
	double C = (r_f.e + (r_z0.e + r_z1.e)) + (r_z2.e + r_z3.e);
	double r_error = (1 + 3 * EPS) * (((1 + 2 * EPS) * (ax / D) - (1 - 3 * EPS) * B) + C);

    return AJRC(r_f.z, r_z0.z, r_z1.z, r_z2.z, r_z3.z, r_error);

}
const AJRC sqrt(const AJRC&x) {

	double xdist = size(x);
	double ax = absUB(x.f);
	double D = ax - (1 + EPS) * (xdist + x.e);
	if(!(D > 0)) {
		return AJRC(0, 0, 0, 0, 0, 0, 0, (1 + 2 * EPS) * sqrt(ax + (xdist + x.e)));
	} else {
		AComplex r_f = sqrt(x.f);
		AComplex t = r_f + r_f;
		AComplex r_z0 = AComplex(x.z0.re, x.z0.im, 0)/t;
		AComplex r_z1 = AComplex(x.z1.re, x.z1.im, 0)/t;
		AComplex r_z2 = AComplex(x.z2.re, x.z2.im, 0)/t;
		AComplex r_z3 = AComplex(x.z3.re, x.z3.im, 0)/t;
	    double C = (r_f.e + (r_z0.e + r_z1.e)) + (r_z2.e + r_z3.e);
		double r_error = (1 + 3 * EPS) * (
		    ((1 + EPS) * sqrt(ax) - (1 - 3 * EPS) * (xdist / (2 * sqrt(ax)) + sqrt(D))) + C);

        return AJRC(r_f.z, r_z0.z, r_z1.z, r_z2.z, r_z3.z, r_error);
	}

}
