#ifndef _Complex_h_
#define _Complex_h_
#include <assert.h>
#include "roundoff.h"

struct XComplex {
	double re;
	double im;
	XComplex(double r = 0,double i = 0):re(r),im(i) {}
};

struct AComplex {
	XComplex z;
	double e;
	AComplex(double r = 0,double i = 0,double err = 0):z(r,i),e(err) {}
};

inline const XComplex operator-(const XComplex& x);
inline const XComplex conj(const XComplex& x);
inline const XComplex re(const XComplex& x);
inline const XComplex im(const XComplex& x);
inline const AComplex operator+(const AComplex& x,const AComplex& y);
inline const AComplex operator+(const XComplex& x,const XComplex& y);
inline const AComplex operator+(const XComplex& x,double y);
inline const AComplex operator-(const AComplex& x,const AComplex& y);
inline const AComplex operator-(const XComplex& x,const XComplex& y);
inline const AComplex operator-(const XComplex& x,double y);
inline const AComplex operator*(const XComplex& x,const XComplex& y);
inline const AComplex operator*(const XComplex& x,double y);
inline const AComplex operator/(const XComplex& x,double y);
inline const double absLB(const XComplex& x);
inline const double absLB(const AComplex& x);
inline const double absUB(const XComplex& x);
inline const double absUB(const AComplex& x);
AComplex operator/(const AComplex& x,const AComplex& y);
AComplex operator/(const XComplex& x,const XComplex& y);
AComplex operator/(double x,const XComplex& y);
AComplex sqrt(const XComplex& x);
#include "Complex.inline"
#endif
