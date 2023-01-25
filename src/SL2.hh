#ifndef _SL2_h_
#define _SL2_h_

#include "types.hh"

template<typename T>
struct SL2 {
  SL2<T>() : a(1), b(0), c(0), d(1) {}
  SL2<T>(
      const T &aa, const T &bb, 
      const T &cc, const T &dd
      ) : a(aa), b(bb), c(cc), d(dd) {}
  T a, b, c, d;
};

template<typename T>
const SL2<T> operator*(const SL2<T> &M, const SL2<T> &N) {
  return SL2<T>(M.a*N.a+M.b*N.c, M.a*N.b+M.b*N.d,
                M.c*N.a+M.d*N.c, M.c*N.b+M.d*N.d);
};

template<typename T>
const SL2<T> inverse(const SL2<T> &M) {
  return SL2<T>(M.d,-M.b,-M.c,M.a);
};

template<typename T>
const T tr(const SL2<T> &M) {
  return M.a + M.d;
};

template<typename T>
const T td(const SL2<T> &M) {
  return M.a - M.d;
};

template<typename T>
const T dist(const SL2<T> &M1, const SL2<T> &M2) {
  return abs(M1.a - M2.a) + abs(M1.b - M2.b) + abs(M1.c - M2.c) + abs(M1.d - M2.d); 
};

template<typename T>
inline const SL2<T> pow(const SL2<T> &M, int n) {
  SL2<T> A; // identity
  if (n == 0) { return A; }
  SL2<T> B;
  if (n < 0) { 
    B = inverse(M);
    n = -n;
  } else {
    B = M;
  }
  while (n > 1) {
    if (n & 1) { // n odd
      A = A*B;
    }
    B = B*B;
    n /= 2; // int division
  } 
  return A*B;
};

template<typename T>
void print_SL2(const SL2<T>& x) {
  print_type("a =", x.a);
  print_type("b =", x.b);
  print_type("c =", x.c);
  print_type("d =", x.d);
}

template<typename T>
const SL2<T> operator-(const SL2<T> &x) {
  return SL2<T>(-x.a, -x.b, -x.c, -x.d);
}

#endif

/* Don't need and not in SL2 
template<typename T>
const SL2<T> operator+(const SL2<T> &x, const SL2<T> &y) {
  return SL2<T>(x.a+y.a, x.b+y.b, x.c+y.c, x.d+y.d);
}


template<typename T>
const SL2<T> operator-(const SL2<T> &x, const SL2<T> &y) {
  return SL2<T>(x.a-y.a, x.b-y.b, x.c-y.c, x.d-y.d);
} */
