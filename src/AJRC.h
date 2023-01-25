#ifndef __AJRC_h_
#define __AJRC_h_
#include "Complex.h"
#include <assert.h>
#include <stdio.h>
#include "roundoff.h"

struct AJRC {
	AJRC(const XComplex& f  = 0,
	   const XComplex& z0 = 0,
	   const XComplex& z1 = 0,
	   const XComplex& z2 = 0,
	   const XComplex& z3 = 0,
	   double err = 0) : f{f},z0{z0},z1{z1},z2{z2},z3{z3},e{err},
		size((1 + 2*EPS) * ((absUB(z0) + absUB(z1)) + (absUB(z2) + absUB(z3)))), 
		size_re((1 + 2*EPS) * 
            ((fabs(z0.re) + fabs(z1.re)) + (fabs(z2.re) + fabs(z3.re)))) {} 
	XComplex f;
	XComplex z0;
	XComplex z1;
	XComplex z2;
	XComplex z3;
	double e;

	double size;
	double size_re;
};

inline const AJRC eye(const AJRC& x) { return AJRC(XComplex(0,1)); };
inline const AJRC operator-(const AJRC& x);
inline const AJRC conj(const AJRC& x);
inline const AJRC re(const AJRC& x);
inline const AJRC im(const AJRC& x);
inline const AJRC operator+(const AJRC& x,const AJRC& y);
inline const AJRC operator-(const AJRC& x,const AJRC& y);
inline const AJRC operator+(const AJRC& x,double y);
inline const AJRC operator-(const AJRC& x,double y);
inline const AJRC operator*(const AJRC& x,double y);
inline const AJRC operator/(const AJRC& x,double y);
inline const AJRC abs_sqrd(const AJRC& x);
inline const AJRC abs(const AJRC& x);
inline const AJRC abs_re(const AJRC& x);
inline const double size(const AJRC& x);
inline const double absUB(const AJRC& x);
inline const double absLB(const AJRC& x);
const AJRC operator*(const AJRC& x,const AJRC& y);
const AJRC operator/(const AJRC& x,const AJRC& y);
const AJRC operator/(double x,const AJRC& y);
const AJRC sqrt(const AJRC& x);

#include "AJRC.inline"
#endif // __AJRC_h
