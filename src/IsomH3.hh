#ifndef __IsomH3_h
#define __IsomH3_h
#include <math.h>
#include <string>
#include "SL2.hh"
#include "Generators.hh"
#include "roundoff.h"
#include "types.hh"
#include "assert.h"

extern bool g_debug;

template<typename T>
inline const T mobius(const SL2<T> &x, const T &p) {
  return ((x.a * p) + x.b) / ((x.c * p) + x.d);
}

template<typename T>
const T four_cosh_re_length(const SL2<T>& w) {
  return abs_sqrd(tr(w)) + abs(tr(w) * tr(w) - 4);
}

template<typename T>
const T cosh_move_j(const SL2<T>& w) {
  T q = abs_sqrd(w.c) + abs_sqrd(w.d);
  T z = w.a * conj(w.c) + w.b * conj(w.d);
  return (abs_sqrd(z) + (q - 1) * (q - 1))/(q * 2) + 1; 
}

template<typename T>
const T two_cosh_dist(T& two_sh_sq_half) {
  return abs(two_sh_sq_half + 2) + abs(two_sh_sq_half);
}

template<typename T>
const T four_cosh_dist(T& four_sh_sq_half) {
  return abs(four_sh_sq_half + 4) + abs(four_sh_sq_half);
}

template<typename T>
const T norm_sqrd(const SL2<T>& w1, const SL2<T>& w2) {
  return (tr(w1) * tr(w1) - 4) * (tr(w2) * tr(w2) - 4);
}

// TODO: Check that we do need the product of the square roots and
// not the square root of the product.
template<typename T>
const T norm(const SL2<T>& w1, const SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T sh1 = sqrt((tr(w1) * tr(w1) - 4));
  T sh2 = sqrt((tr(w2) * tr(w2) - 4));
  // TODO this might not be reliable
  if (absUB(tr1 + sh1) < 2) {
    if (g_debug) {
      fprintf(stderr, "Flipping sqrt sign (norm): ");
      T t = tr1 + sh1;
      print_type(t);
    }
    sh1 = -sh1;
  }
  // TODO this might not be reliable
  if (absUB(tr2 + sh2) < 2) {
    if (g_debug) {
      fprintf(stderr, "Flipping sqrt sign (norm): ");
      T t = tr2 + sh2;
      print_type(t);
    }
    sh2 = -sh2;
  }
  return sh1 * sh2;
}

template<typename T>
const T cosh_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  return td(w1) * td(w2) + (w1.b * w2.c + w1.c * w2.b) * 2;
}

template<typename T>
const T cosh_perp(const SL2<T>& w1, const SL2<T>& w2) {
  return cosh_perp_normed(w1,w2)/norm(w1,w2);
}

template<typename T>
const T cosh_perp_sqrd_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T z = cosh_perp_normed(w1,w2); 
  return z*z; 
}

template<typename T>
const T cosh_perp_sqrd(const SL2<T>& w1, const SL2<T>& w2) {
  return cosh_perp_sqrd_normed(w1,w2)/norm_sqrd(w1,w2);
}

template<typename T>
const T abs_cosh_perp_sqrd_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T z = cosh_perp_normed(w1,w2); 
  return abs_sqrd(z); 
}

template<typename T>
const T abs_cosh_perp_sqrd(const SL2<T>& w1, const SL2<T>& w2) {
  // TODO: check if this gives best error as 
  // both sqrt and division make error blow up
  return abs_cosh_perp_sqrd_normed(w1,w2)/abs(norm_sqrd(w1,w2));
}

template<typename T>
const T sinh_perp_sqrd_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T ch = cosh_perp_normed(w1,w2); 
  return ch * ch - norm_sqrd(w1,w2);
}

template<typename T>
const T sinh_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T ch = cosh_perp_normed(w1,w2); 
  T n_sqrd = norm_sqrd(w1,w2);
  T sh = sqrt(ch * ch - n_sqrd);
  // TODO this might not be reliable
  if (absUB(ch + sh) < absUB(norm(w1,w2))) {
    if (g_debug) {
      fprintf(stderr, "Flipping sqrt sign (sinh_perp_normed): ");
      T t = ch + sh;
      print_type(t);
    }
    sh = -sh;
  }
  return sh;
}

template<typename T>
const T sinh_perp(const SL2<T>& w1, const SL2<T>& w2) {
  T ch = cosh_perp(w1,w2);
  T sh = sqrt(ch * ch - 1);
  // TODO this might not be reliable
  if (absUB(ch + sh) < 1) {
    if (g_debug) {
      fprintf(stderr, "Flipping sqrt sign (sinh_perp): ");
      T t = ch + sh;
      print_type(t);
    }
    sh = -sh;
  }
  return sh;
}

template<typename T>
const T abs_sinh_perp_sqrd_normed(const SL2<T>& w1, const SL2<T>& w2) {
  return abs(sinh_perp_sqrd_normed(w1,w2));
}

template<typename T>
const T abs_sinh_perp_sqrd(const SL2<T>& w1, const SL2<T>& w2) {
  // TODO: check if this gives best error as 
  // both sqrt and division make error blow up
  return abs_sinh_perp_sqrd_normed(w1,w2)/abs(norm_sqrd(w1,w2));
}

template<typename T>
const T cosh_2_re_perp(const SL2<T>& w1, const SL2<T>& w2) {
  T abs_cp_sqrd = abs_cosh_perp_sqrd(w1,w2);
  T abs_sp_sqrd = abs_sinh_perp_sqrd(w1,w2);
  return abs_cp_sqrd + abs_sp_sqrd;
}

template<typename T>
const T cosh_2_re_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T abs_cp_sqrd_normed = abs_cosh_perp_sqrd_normed(w1,w2);
  T abs_sp_sqrd_normed = abs_sinh_perp_sqrd_normed(w1,w2);
  return abs_cp_sqrd_normed + abs_sp_sqrd_normed;
}

template<typename T>
const T sinh_2_re_perp(const SL2<T>& w1, const SL2<T>& w2) {
  T ch = cosh_2_re_perp(w1,w2);
  T sh = sqrt(ch * ch - 1);
  // TODO this might not be reliable
  if (absUB(ch + sh) < 1) {
    if (g_debug) {
      fprintf(stderr, "Flipping sqrt sign (sinh_2_re_perp): ");
      T t = ch + sh;
      print_type(t);
    } 
    sh = -sh;
  }
  return sh;
}

template<typename T>
const T sinh_2_re_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T ch = cosh_2_re_perp_normed(w1,w2);
  T n = norm_sqrd(w1,w2);
  T abs_n_sqrd = abs_sqrd(n);
  T sh = sqrt(ch * ch - abs_n_sqrd);
  // TODO this might not be reliable
  if (absUB(ch + sh) < absUB(n)) {
    if (g_debug) {
      fprintf(stderr, "Flipping sqrt sign (sinh_2_re_perp_normed): ");
      T t = ch + sh;
      print_type(t);
    } 
    sh = -sh;
  }
  return sh;
}

template<typename T>
const T e_perp(const SL2<T>& w1, const SL2<T>& w2) {
  T cp = cosh_perp(w1,w2);
  T sp = sinh_perp(w1,w2);
  return cp + sp;
}

template<typename T>
const T e_2_re_perp(const SL2<T>& w1, const SL2<T>& w2) {
  return abs_sqrd(e_perp(w1,w2));
}

template<typename T>
const T e_2_re_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T cp = cosh_perp_normed(w1,w2);
  T sp = sinh_perp_normed(w1,w2);
  return abs_sqrd(cp + sp);
}

template<typename T>
const T e_m_perp(const SL2<T>& w1, const SL2<T>& w2) {
  T cp = cosh_perp(w1,w2);
  T sp = sinh_perp(w1,w2);
  return cp - sp;
}

template<typename T>
const T e_m_2_re_perp_normed(const SL2<T>& w1, const SL2<T>& w2) {
  T cp = cosh_perp_normed(w1,w2);
  T sp = sinh_perp_normed(w1,w2);
  return abs_sqrd(cp - sp);
}

template<typename T>
const double e_re_perp_UB(const SL2<T>& w1, const SL2<T>& w2) {
  return absUB(e_perp(w1,w2));
}

template<typename T>
const double e_re_perp_LB(const SL2<T>& w1, const SL2<T>& w2) {
  return absLB(e_perp(w1,w2));
}

template<typename T>
const std::pair<T,T> four_cosh_margulis_simple(
    const SL2<T>& w1, const SL2<T>& w2) {
  /*
     print_center("w1.a:", w1.a);
     print_center("w1.b:", w1.b);
     print_center("w1.c:", w1.c);
     print_center("w1.d:", w1.d);
     print_center("w2.a:", w2.a);
     print_center("w2.b:", w2.b);
     print_center("w2.c:", w2.c);
     print_center("w2.d:", w2.d);
   */
  // retuns 4 cosh( margulis ) and exp(2t) for w1,w2
  // TODO: check that the words don't commute
  const T tr1 = w1.a + w1.d;
  const T tr2 = w2.a + w2.d;
  const T x1 = abs_sqrd(tr1);
  const T x2 = abs(tr1*tr1 - 4);
  const T x2_sqrd = abs_sqrd(tr1*tr1 - 4);
  const T y1 = abs_sqrd(tr2);
  const T y2 = abs(tr2*tr2 - 4);
  const T y2_sqrd = abs_sqrd(tr2*tr2 - 4);
  // exp(2re(P))x2y2
  const T e_2_re_perp_n = e_2_re_perp_normed(w1,w2);
  const T e_m_2_re_perp_n = e_m_2_re_perp_normed(w1,w2);
  const T e_2_re_p = e_2_re_perp(w1,w2);
  // cosh(2(re(P))x2y2
  const T ch_2_re_perp_n = cosh_2_re_perp_normed(w1,w2); 
  // sinh(2(re(P))x2y2
  const T sh_2_re_perp_n = sinh_2_re_perp_normed(w1,w2); 

  /*
     printf("***********************************\n");
     print_center("tr(w1):", tr1);
     print_center("tr(w2):", tr2);
     print_center("x1:", x1);
     print_center("x2:", x2);
     print_center("y1:", y1);
     print_center("y2:", y2);
     print_center("cosh(perp)sqrt((tr1^2-4)(tr2^2-4)):",
        cosh_perp_normed(w1,w2));
     print_center("sinh(perp)sqrt((tr1^2-4)(tr2^2-4)):",
        sinh_perp_normed(w1,w2));
     print_center("exp(2re(perp))x2y2:", e_2_re_perp_n);
     print_center("exp(-2re(perp))x2y2:", e_m_2_re_perp_n);
     print_center("cosh(2re(perp))x2y2:", ch_2_re_perp_n);
     print_center("sinh(2re(perp))x2y2:", sh_2_re_perp_n);
     printf("***********************************\n");
     print_center("e_2_re_perp_n*y1 - x1*y2_sqrd:",
        e_2_re_perp_n*y1 - x1*y2_sqrd);
     print_center("e_2_re_perp_n - y2_sqrd:", e_2_re_perp_n - y2_sqrd);
     print_center("a/b:",
        (e_2_re_perp_n*y1 - x1*y2_sqrd)/(e_2_re_perp_n - y2_sqrd));
     print_center("sinh(2re(perp))x2y2:", sh_2_re_perp_n);
     print_center("(y1-x1)+sqrt((y1-x1)*(y1-x1) +
        ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd)):",
        (y1-x1)+sqrt((y1-x1)*(y1-x1) +
        ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd)));
     print_center("a/b:",
        sh_2_re_perp_n/((y1-x1)+sqrt((y1-x1)*(y1-x1) +
        ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd))));
     print_center("x2*((y1-x1) + sqrt((y1-x1)*(y1-x1) +
        ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd))):",
        x2*((y1-x1) + sqrt((y1-x1)*(y1-x1) +
        ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd))));
     print_center("x2_sqrd - e_m_2_re_perp_n:", x2_sqrd - e_m_2_re_perp_n);
     print_center("a/b:",
        (x2*((y1-x1) + sqrt((y1-x1)*(y1-x1) +
        ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd)))) /
        (x2_sqrd - e_m_2_re_perp_n));
     printf("***********************************\n");
   */

  std::vector<T> versions;
  versions.push_back((e_2_re_perp_n*y1 - x1*y2_sqrd) / 
      (e_2_re_perp_n - y2_sqrd) +
      sh_2_re_perp_n/((y1-x1)+sqrt((y1-x1)*(y1-x1) +
          ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd))));

  // swap x and y to hope for better error
  // (formula is symmetric, but not in an obious way)
  versions.push_back((e_2_re_perp_n*x1 - y1*x2_sqrd) /
      (e_2_re_perp_n - x2_sqrd) +
      sh_2_re_perp_n/((x1-y1) + sqrt((x1-y1)*(x1-y1) +
          ((ch_2_re_perp_n*2 - y2_sqrd) - x2_sqrd))));

  // Alternate version where first denominator
  // is scaled by e^(-2re(P))y2/x2.
  // The current version should keep the denominator larger
  //  const T four_cosh_marg_v1_alt1 = (y1*x2_sqrd - x1*e_m_2_re_perp_n) /
  //  (x2_sqrd - e_m_2_re_perp_n) +
  //             sh_2_re_perp_n/((y1-x1)+sqrt((y1-x1)*(y1-x1) +
  //             ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd)));

  // Alternate version where we split off x1.
  // The current version should keep error slightly smaller
  // (since fewer additions)
  //  const T four_cosh_marg_v1_alt2 = (x1 + ((y1-x1)*x2_sqrd) /
  //  (x2_sqrd - e_m_2_re_perp_n)) +
  //             sh_2_re_perp_n/((y1-x1)+sqrt((y1-x1)*(y1-x1) +
  //             ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd)));

  // Collect everything as one fraction.
  // Helps when re(P) is very close to zero 
  versions.push_back(x1 + ((y1-x1)*(e_2_re_perp_n-x2_sqrd) -
        sh_2_re_perp_n*((y1-x1) -
          sqrt((y1-x1)*(y1-x1) +
            ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd)))) /
      ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd));

  // swap x and y to hope for better error
  // (formula is symmetric, but not in an obious way)
  versions.push_back(y1 + ((x1-y1)*(e_2_re_perp_n-y2_sqrd)
        - sh_2_re_perp_n*((x1-y1) -
          sqrt((x1-y1)*(x1-y1) +
            ((ch_2_re_perp_n*2 - y2_sqrd) - x2_sqrd)))) /
      ((ch_2_re_perp_n*2 - y2_sqrd) - x2_sqrd));

  /* 
     print_type("4cosh(margulis) v1:", versions[0]);
     print_type("4cosh(margulis) v2:", versions[1]);
     print_type("4cosh(margulis) v3:", versions[2]);
     print_type("4cosh(margulis) v4:", versions[3]);
   */

  std::sort(versions.begin(), versions.end(), sort_comp<T>);

  T four_cosh_marg = versions[0];

  T exp_2_t = (x2*((y1-x1) + sqrt((y1-x1)*(y1-x1) +
          ((ch_2_re_perp_n*2 - x2_sqrd) - y2_sqrd)))) /
    (x2_sqrd - e_m_2_re_perp_n);
  /*
     print_center("exp(2re(perp)):", e_2_re_p);
     print_center("exp(2t):", (const T) exp_2_t);
     print_center("exp(2t) - exp(2re(perp)):", exp_2_t - e_2_re_p);
     print_center("exp(2re(perp)) exp(2t):", e_2_re_p - exp_2_t);
   */

  // if we realize the margulis number outside the ortho segement,
  // it means it is on one of the end point
  const T one = T(1);
  const T zero = T(0);
  T fc_marg_lb = max_local(x2 + x1, y2 + y1);
  // print_type("lower bound for four cosh margulis:", fc_marg_lb);
  if (absUB(exp_2_t) < 1.0) {
    four_cosh_marg = x2 + x1;
    exp_2_t = one;
  } else if (absLB(exp_2_t - 1.0) == 0) {
    four_cosh_marg = max_local(four_cosh_marg, fc_marg_lb);
    exp_2_t = max_local(exp_2_t, one);
  } else if (absLB(exp_2_t) > absUB(e_2_re_p)) {
    four_cosh_marg = y2 + y1;
    exp_2_t = e_2_re_p;  
  } else if (absLB(exp_2_t - e_2_re_p) == 0) {
    four_cosh_marg = max_local(four_cosh_marg, fc_marg_lb);
    exp_2_t = min_local(exp_2_t, e_2_re_p);
  }
  if (absLB(y2 + y1) > absUB((e_2_re_p*x2)+x1)) { 
    four_cosh_marg = y2+y1;
    exp_2_t = e_2_re_p; 
  }
  if (absLB(x2 + x1) > absUB((e_2_re_p*y2)+y1)) { 
    four_cosh_marg = x2+x1;
    exp_2_t = zero; 
  }

  //  print_type("4cosh(margulis) final:", four_cosh_marg);
  //  print_type("exp_2_t final:", exp_2_t);

  std::pair<T,T> result(four_cosh_marg, exp_2_t);

  return result;
}

// We compute |tr(w1)^2 - 4| + |tr(w1 w2 W1 W2) - 2|
// with optimzation for x and y specifically
template<typename T>
const T jorgensen(const SL2<T>& w1, const SL2<T>& w2) {
  SL2<T> W1 = inverse(w1); 
  SL2<T> W2 = inverse(w2);
  SL2<T> comm = w1*w2*W1*W2;
  return abs(tr(w1) * tr(w1) - 4) + abs(tr(comm) - 2);
}

template<typename T>
const std::pair<T,T> fixed_points(const SL2<T>& w) {
  std::pair<T,T> result(
      (td(w) + sqrt(tr(w) * tr(w) - 4)) / (w.c * 2),
      (td(w) - sqrt(tr(w) * tr(w) - 4)) / (w.c * 2));
  return result;
}

// Complex distance between {zm, zp} and {0, infty} 
template<typename T>
const T sinh_sqrd_half_perp_zero_inf(T& zm, T& zp) {
  return zm / (zp - zm); 
}

// Complex distance between {zm, zp} and {-1, 1} 
template<typename T>
const T sinh_sqrd_half_perp_mp_one(T& zm, T& zp) {
  T one = T(1);
  return (((zm + one) * (zp - one)) * 0.5) / (zm - zp) ; 
}

// Complex distance between {zm, zp} and {-1, 1} 
template<typename T>
const T sinh_sqrd_half_perp_mp_iye(T& zm, T& zp) {
  T one = T(1);
  T iye = eye(one);
  return (((one - iye * zm) * (zp - iye)) * 0.5) / (zm - zp) ; 
}

#endif // __IsomH3_h

