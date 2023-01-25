#ifndef __AJCC_h_
#define __AJCC_h_
#include "Complex.h"
#include <assert.h>
#include <stdio.h>
#include "roundoff.h"

// note: costs (1+EPS/2) that must be computed later
inline double abs_taxi(const XComplex& x) {
  return fabs(x.re) + fabs(x.im);
}

struct AJCC {
	AJCC(const XComplex& f  = 0,
	   const XComplex& z0 = 0,
	   const XComplex& z1 = 0,
	   const XComplex& w0 = 0,
	   const XComplex& w1 = 0,
	   double err = 0) : f{f},z0{z0},z1{z1},w0{w0},w1{w1},e{err},
		size((1+3*EPS) * ((abs_taxi(z0) + abs_taxi(z1)) +
                          (abs_taxi(w0) + abs_taxi(w1)))) {} 
	XComplex f;
	XComplex z0;
	XComplex z1;
	XComplex w0;
	XComplex w1;
	double e;

	double size;
};

inline const AJCC eye(const AJCC& x) { return AJCC(XComplex(0,1)); };
inline const AJCC operator-(const AJCC& x);
inline const AJCC conj(const AJCC& x);
inline const AJCC re(const AJCC& x);
inline const AJCC im(const AJCC& x);
inline const AJCC operator+(const AJCC& x,const AJCC& y);
inline const AJCC operator-(const AJCC& x,const AJCC& y);
inline const AJCC operator+(const AJCC& x,double y);
inline const AJCC operator-(const AJCC& x,double y);
inline const AJCC operator*(const AJCC& x,double y);
inline const AJCC operator/(const AJCC& x,double y);
inline const double absUB(const AJCC& x);
inline const double absLB(const AJCC& x);
inline const AJCC abs_sqrd(const AJCC& x);
inline const AJCC abs(const AJCC& x);
inline const double size(const AJCC& x);
const AJCC operator*(const AJCC& x,const AJCC& y);
const AJCC operator/(const AJCC& x,const AJCC& y);
const AJCC operator/(double x,const AJCC& y);
const AJCC sqrt(const AJCC& x);

#include "AJCC.inline"
#endif // __AJCC_h
