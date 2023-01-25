#ifndef __Params_h
#define __Params_h
#include <math.h>
#include <string>
#include "SL2.hh"
#include "Generators.hh"
#include "roundoff.h"
#include "types.hh"
#include "assert.h"

// Eliminate bad boxes that can't gerate non-elementay groups
template<typename T>
const T jorgensen_xy(const Params<T>& p) {
  T z = p.sinhLy2 * p.sinhperp; 
  return (abs_sqrd(z) + 1) * abs_sqrd(p.sinhLx2) * 4;
}

// Eliminate bad boxes that can't gerate non-elementay groups
template<typename T>
const T jorgensen_yx(const Params<T>& p) {
  T z = p.sinhLx2 * p.sinhperp; 
  return (abs_sqrd(z) + 1) * abs_sqrd(p.sinhLy2) * 4;
}

// Complex distance between axis(g) and w(axis(g)) 
template<typename T>
const T four_sinh_sqrd_half_perp(const SL2<T>& w,
    const Params<T>& p, const char g) {
  T b_coeff, c_coeff, shL2;
  if (g == 'x') {
    shL2 = p.sinhLx2;
    b_coeff = p.expdx;
    c_coeff = p.expmdx;
  } else {
    shL2 = p.sinhLy2;
    b_coeff = p.expmdyf;
    c_coeff = p.expdyf;
 }
  T z = w.c * c_coeff - w.b * b_coeff;
  // formula by using crossratios
  return td(w) * td(w) - z * z;  
}

template<typename T>
const T jorgensen_gw(const SL2<T>& w, const Params<T>& p, const char g) {
  T shL2;
  if (g == 'x') {
    shL2 = p.sinhLx2;
  } else {
    shL2 = p.sinhLy2;
 }
  return (abs(four_sinh_sqrd_half_perp(w, p, g)) + 4) * abs_sqrd(shL2);
}

template<typename T>
const T jorgensen_wg(const SL2<T>& w, const Params<T>& p, const char g) {
  T shL2;
  if (g == 'x') {
    shL2 = p.sinhLx2;
  } else {
    shL2 = p.sinhLy2;
 }
  return abs(tr(w) * tr(w) - 4) + 
    abs(four_sinh_sqrd_half_perp(w, p, g)) * abs_sqrd(shL2);
}

// Distance between axis(x) and w(axis(x)) 
template<typename T>
const T four_cosh_dist(const SL2<T>& w, const Params<T>& p, const char g) {
  T four_sh_sq_half_perp = four_sinh_sqrd_half_perp(w, p, g); 
  return four_cosh_dist(four_sh_sq_half_perp);
}

// Complex distance between axis(x) and w(axis(y)) 
template<typename T>
const T four_sinh_sqrd_half_perp_xwy(const SL2<T>& w, const Params<T>& p) {
  T z = ((w.a * w.a) * p.expdyf - (w.b * w.b) * p.expmdyf) * p.expdx
      + ((w.d * w.d) * p.expmdyf - (w.c * w.c) * p.expdyf) * p.expmdx;
  return z - 2;
}

// Distance between axis(x) and w(axis(y)) 
template<typename T>
const T four_cosh_dist_xwy(const SL2<T>& w, const Params<T>& p) {
  T z = ((w.a * w.a) * p.expdyf - (w.b * w.b) * p.expmdyf) * p.expdx
      + ((w.d * w.d) * p.expmdyf - (w.c * w.c) * p.expdyf) * p.expmdx;
  return  abs(z - 2) + abs(z + 2);
}

// Complex distance between w(axis(g)) and {0, infty} 
template<typename T>
const T two_sinh_sqrd_half_perp_wg_zero_inf(const SL2<T>& w,
    const Params<T>& p, const char g) {
  T ac_coeff, bd_coeff;
  if (g == 'x') {
    ac_coeff = p.expmdx;
    bd_coeff = p.expdx;
  } else {
    ac_coeff = p.expdyf;
    bd_coeff = p.expmdyf;
 }
  T one = T(1);
  return (w.b * w.d) * bd_coeff - (w.a * w.c) * ac_coeff - one; 
}

template<typename T>
const T cosh_sqrd_perp_wg_zero_inf(const SL2<T>& w,
    const Params<T>& p, const char g) {
  T ac_coeff, bd_coeff;
  if (g == 'x') {
    ac_coeff = p.expmdx;
    bd_coeff = p.expdx;
  } else {
    ac_coeff = p.expdyf;
    bd_coeff = p.expmdyf;
 }
  T z = (w.b * w.d) * bd_coeff - (w.a * w.c) * ac_coeff;
  return z * z; 
}

template<typename T>
const T sinh_sqrd_perp_wg_zero_inf(const SL2<T>& w,
    const Params<T>& p, const char g) {
  T ac_coeff, bd_coeff;
  if (g == 'x') {
    ac_coeff = p.expmdx;
    bd_coeff = p.expdx;
  } else {
    ac_coeff = p.expdyf;
    bd_coeff = p.expmdyf;
 }
  T ad = w.a * w.d;
  T bc = w.b * w.c;
  T aco = (w.a * w.c) * ac_coeff;
  T bdo = (w.b * w.d) * bd_coeff;
  return (aco * aco + bdo * bdo) - (ad * ad + bc * bc); 
}

// Complex distance between w(axis(g)) and {-1, 1} 
template<typename T>
const T four_sinh_sqrd_half_perp_wg_mp_one(const SL2<T>& w,
    const Params<T>& p, const char g) {
  T ac_coeff, bd_coeff;
  if (g == 'x') {
    ac_coeff = p.expmdx;
    bd_coeff = p.expdx;
  } else {
    ac_coeff = p.expdyf;
    bd_coeff = p.expmdyf;
 }
  T two = T(2);
  return (w.d * w.d - w.b * w.b) * bd_coeff
    - (w.c * w.c - w.a * w.a) * ac_coeff
    - two;
}

// Complex distance between w(axis(x)) and {-I, I} 
template<typename T>
const T four_sinh_sqrd_half_perp_wg_mp_iye(const SL2<T>& w,
    const Params<T>& p, const char g) {
  T iye = eye(T(1));
  T ac_coeff, bd_coeff;
  if (g == 'x') {
    ac_coeff = p.expmdx;
    bd_coeff = p.expdx;
  } else {
    ac_coeff = p.expdyf;
    bd_coeff = p.expmdyf;
 }
  T two = T(2);
  return ((w.d * w.d + w.b * w.b) * bd_coeff
    - (w.c * w.c + w.a * w.a) * ac_coeff) * iye
    - two;
}

#endif // __Params_h
