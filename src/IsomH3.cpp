#include "IsomH3.hh"
#include "roundoff.h"

// TODO: Check that we do need the product of the square roots and
// note the square root of the product.
template<typename T>
const T norm(const SL2<T>& w1, SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  return sqrt(tr1*tr1-4)*sqrt(tr2*tr2-4);
}

template<typename T>
const T norm_sqrd(const SL2<T>& w1, SL2<T>& w2) {
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  return (tr1*tr1-4)*(tr2*tr2-4);
}

template<typename T>
const T cosh_perp_normed(const SL2<T>& w1, SL2<T>& w2) {
  T td1 = w1.a - w1.d;
  T td2 = w2.a - w2.d;
  return td1 * td2 + 2*(w1.b * w2.c + w1.c * w2.b);
}

template<typename T>
const T cosh_perp(const SL2<T>& w1, SL2<T>& w2) {
  T z = cosh_perp_normed(w1,w2); 
  return z/norm(w1,w2);
}

template<typename T>
const T cosh_perp_sqrd_normed(const SL2<T>& w1, SL2<T>& w2) {
  T z = cosh_perp_normed(w1,w2); 
  return z*z; 
}

template<typename T>
const T cosh_perp_sqrd(const SL2<T>& w1, SL2<T>& w2) {
  return cosh_perp_sqrd_normed(w1,w2)/norm_sqrd(w1,w2);
}

template<typename T>
const T abs_cosh_perp_sqrd_normed(const SL2<T>& w1, SL2<T>& w2) {
  T z = cosh_perp_normed(w1,w2); 
  return abs_sqrd(z); 
}

template<typename T>
const T abs_cosh_perp_sqrd(const SL2<T>& w1, SL2<T>& w2) {
  return abs_cosh_perp_sqrd_normed(w1,w2)/abs(norm_sqrd(w1,w2)); // TODO: check if this gives best error as both sqrt and division make error blow up
}

template<typename T>
const T sinh_perp_sqrd_normed(const SL2<T>& w1, SL2<T>& w2) {
  T z = cosh_perp_normed(w1,w2); 
  return z*z - norm_sqrd(w1,w2);
}

//template<typename T>
//const T sinh_perp_normed(const SL2<T>& w1, SL2<T>& w2) {
//  return sqrt(cosh_perp_sqrd(w1,w2)-1)*norm(w1,w2);
//}

template<typename T>
const T sinh_perp(const SL2<T>& w1, SL2<T>& w2) {
  T z = cosh_perp(w1,w2); 
  return sqrt(z*z-1); // TODO: check if this gives best error
}

template<typename T>
const T abs_sinh_perp_sqrd_normed(const SL2<T>& w1, SL2<T>& w2) {
  return abs(sinh_perp_sqrd_normed(w1,w2));
}

template<typename T>
const T abs_sinh_perp_sqrd(const SL2<T>& w1, SL2<T>& w2) {
  return abs_sinh_perp_sqrd_normed(w1,w2)/abs(norm_sqrd(w1,w2)); // TODO: check if this gives best error
}

template<typename T>
const T cosh_2_re_perp(const SL2<T>& w1, SL2<T>& w2) {
  T abs_cp_sqrd = abs_cosh_perp_sqrd(w1,w2);
  T abs_sp_sqrd = abs_sinh_perp_sqrd(w1,w2);
  return abs_cp_sqrd + abs_sp_sqrd;
}

template<typename T>
const T sinh_2_re_perp(const SL2<T>& w1, SL2<T>& w2) {
  T c2rep = cosh_2_re_perp(w1,w2);
  return sqrt(c2rep*c2rep - 1);
}

template<typename T>
const double four_cosh_margulis(const SL2<T>& w1, const SL2<T>& w2, bool upper) {
  // retuns 4 cosh( margulis ) for w1,w2
  // TODO : Check signs ERROR ROUNDING
  T tr1 = w1.a + w1.d;
  T tr2 = w2.a + w2.d;
  T x1 = tr1*tr1;
  T x2 = tr1*tr1 - 4;
  T y1 = tr2*tr2;
  T y2 = tr2*tr2 - 4;
  // Normed cosh and sihn values
  T cp = cosh_perp(w1,w2);
  T sp = sinh_perp(w1,w2);
  T ep = cp + sp;
  // exp(2re(P)) 
  T e_2_re_perp = abs_sqrd(ep);
  // cosh(2(re(P))
  T ch_2_re_perp = cosh_2_re_perp(w1,w2); 
  // sinh(2(re(P))
  T sh_2_re_perp = sinh_2_re_perp(w1,w2); 

  T result = (e_2_re_perp*(y1*x2) - (x1*y2))/(e_2_re_perp*x2 - y2) +
             (sh_2_re_perp*(x2*y2))/((y1-x1)+sqrt((y1-x1)^2 + (ch_2_re_perp*2)*(x2*y2) - (x2*x2) - (y2*y2)));

  if (upper) {
    return absUB(result);
  } else {
    return absLB(result);
  }
} 

template<typename T>
const T cosh_2_perp(const SL2<T>& w1, SL2<T>& w2) {
  return 2 * cosh_perp_sqrd(w1,w2) - 1;
}

template<typename T>
const T cosh_perp_x(const SL2<T>& w) {
  // returns cosh(R) where R is the complex distance from axis(w) to (0,inf)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  return td/sqrt(tr*tr - 4); 
}

template<typename T>
const T cosh_perp_x_normed(const SL2<T>& w, Params<T>& params) {
  T td = w.a - w.d;
  return 2*td*params.sinhL2; 
}

template<typename T>
const T sinh_perp_x_normed(const SL2<T>& w, Params<T>& params) {
  return 4*sqrt(1 - w.a * w.d)*params.sinhL2; 
}

template<typename T>
const T cosh_perp_x_sqrd(const SL2<T>& w) {
  // return cosh(R)^2 where R is the complex distance from axis(w) to (0,inf)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  return (td*td)/(tr*tr - 4); 
}

template<typename T>
const T cosh_perp_x_sqrd_normed(const SL2<T>& w, Params<T>& params) {
  T td = w.a - w.d;
  return 4*(td*td)*(params.coshL2*params.coshL2-1); 
}

template<typename T>
const T sinh_perp_x_sqrd_normed(const SL2<T>& w, Params<T>& params) {
  return 16*(1 - w.a * w.d)*(params.coshL2*params.coshL2-1); 
}

template<typename T>
const T cosh_perp_y(const SL2<T>& w, Params<T>& params) {
  // returns cosh(R) where R is the complex distance from axis(w) to axis(y)
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  return (td * params.coshP + dd * params.sinhP)/sqrt(tr*tr - 4); 
}

template<typename T>
const T cosh_perp_y_normed(const SL2<T>& w, Params<T>& params) {
  T td = w.a - w.d;
  T dd = w.b - w.c;
  return 2*(td * params.coshP + dd * params.sinhP)*params.sinhD2;
}

template<typename T>
const T sinh_perp_y_normed(const SL2<T>& w, Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T z  = td * params.coshP + dd * params.sinhP;
  return 2*sqrt(z*z-(tr*tr - 4))*params.sinhD2;
}

template<typename T>
const T cosh_perp_y_sqrd(const SL2<T>& w, Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T z  = td * params.coshP + dd * params.sinhP;
  return (z*z)/(tr*tr - 4); 
}

template<typename T>
const T cosh_perp_y_sqrd_normed(const SL2<T>& w, Params<T>& params) {
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T z  = td * params.coshP + dd * params.sinhP;
  return 4*(z*z)*(params.coshD2-1);
}

template<typename T>
const T sinh_perp_y_sqrd_normed(const SL2<T>& w, Params<T>& params) {
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T dd = w.b - w.c;
  T z  = td * params.coshP + dd * params.sinhP;
  return 4*((z*z)-(tr*tr - 4))*(params.coshD2*params.coshD2-1);
}

//template<typename T>
//const double cosh_2_re_tube_UB(const SL2<T>& w1, const SL2<T>& w2) {
//} 

/*
const ACJpair fixed_points(const SL2ACJ& w) {
  // Returns the pair of fixed point of the Mobius transform w
  ACJ tr = w.a + w.d;
  ACJ td = w.a - w.d;
  return ACJpair((td + sqrt(tr*tr - 4))/(2*w.c),
                 (td - sqrt(tr*tr - 4))/(2*w.c));
}

const ACJ cosh_perp(const ACJpair& p1, ACJpair& p2) {
  return  ((p1.x + p1.y)(p2.x + p2.y) - 2 * p1.x * p2.y - 2 * p2.x * p2.y)/((p1.x - p1.y)(p2.x - p2.y));
}

const std::Cpair fixed_points(const SL2C& w); 
  // Returns the pair of fixed point of the Mobius transform w
  XComplex tr = (w.a + w.d).z;
  XComplex s = (w.a - w.d).z;
  return Cpair(((s + sqrt(((tr*tr)).z-4).z)/(2*w.c)).z,
               ((s - sqrt(((tr*tr)).z-4).z)/(2*w.c)).z);
}*/

