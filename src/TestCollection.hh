#ifndef _TestCollection_ 
#define _TestCollection_
#include <unordered_map>
#include <string>
#include <vector>
#include "types.hh"
#include "Box.h"
#include "SL2.hh"
#include "IsomH3.hh"
#include "Params.hh"
#include "RelatorTest.hh"

extern bool g_debug;

extern double g_cosh_marg_upper;
extern double g_cosh_marg_lower;
extern double g_sinh_r;
extern double g_cosh_r;

extern double g_cosh_sym_marg;
extern double g_sinh_sym_r;
extern double g_cosh_sym_r;
extern double g_cosh_sym_2r;

struct RelatorTest;

struct TestCollection {
  int size();
  box_state evaluate_center(int index, Box& box);
  TestResult evaluate_box(int index, Box& box);
  TestResult evaluate_AJ(word_pair& pair, Box& box);
  TestResult evaluate_qrs(Box& box);
  const std::string get_name(int index);
  word_pair get_pair(int index);
  int add(word_pair pair);
  int add(std::string pair);
  void load(const char* file_path);
  void load_relator_test(const char* impos_path,
                         const char* bad_rel_path);
  RelatorTest *relator_test;
  std::map<std::string, int> seen_words;
private:
  word_pair parse_word_pair(std::string buf);
  std::map<word_pair, int> pair_index;
  std::vector<word_pair> pair_vector;
  box_state evaluate_approx(word_pair pair, const Box& params);
  bool ready_for_elliptics_test(SL2<AJ>& w);
  bool only_elliptics(SL2<AJ>& w, Params<AJ>& params);
};

template<typename T>
inline const bool inside_var_nbd(const SL2<T>& w1, const SL2<T>& w2) {
  // Show that w1 and w2 are commuting non-parabolics when discrete
  // So either we have a relator or they are elliptic
  return (absUB(jorgensen(w1, w2)) < 1 || absUB(jorgensen(w2, w1)) < 1)
    && (not_parabolic(w1) || not_parabolic(w2));
}

template<typename T>
inline const bool inside_var_nbd_ne(const SL2<T>& w1, const SL2<T>& w2) {
  // Show that w1 and w2 are commuting loxodromics when discrete,
  // so we have a relator
  // Note,we must test both as elliptic can commute with loxodromic
  return (absUB(jorgensen(w1, w2)) < 1 || absUB(jorgensen(w2, w1)) < 1) 
    && (not_elliptic_or_parabolic(w1) && not_elliptic_or_parabolic(w2));
}

template<typename T>
inline const bool not_parabolic(const SL2<T>& w) {
  return absLB(im(tr(w))) > 0 || (absLB(re(tr(w)) - 2) > 0 
      && absLB(re(tr(w)) + 2) > 0);
}

template<typename T>
inline const bool not_elliptic_or_parabolic(const SL2<T>& w) {
  return absLB(im(tr(w))) > 0 || absLB(re(tr(w))) > 2; 
}

template<typename T>
inline const bool not_identity(const SL2<T>& w) {
  return absLB(w.b) > 0 ||  absLB(w.c) > 0 ||
    ((absLB(w.a - 1) > 0 || absLB(w.d - 1) > 0) 
     && (absLB(w.a + 1) > 0 || absLB(w.d + 1) > 0));
}

template<typename T>
inline const bool tube_hits_axis_two(const T& two_sh_sq_hf_p,
    const T& two_ch_re_tube) {
  T t_ch_d = two_cosh_dist(two_sh_sq_hf_p);
  // sinh(I Pi/4)^2 = -1/2 which means axes meet orthogonally
  return strictly_pos(two_ch_re_tube - t_ch_d)
    && absLB(two_sh_sq_hf_p + 1) > 0;
}

template<typename T>
inline const bool wg_hits_sym_axis(const SL2<T>& w,
    const Params<T>& p, const char g) {
  T t_sh_sq_hf_p = two_sinh_sqrd_half_perp_wg_zero_inf(w, p, g);
  return tube_hits_axis_two(t_sh_sq_hf_p, p.coshdx * 2);
}

template<typename T>
inline const bool w_in_sym_search(const SL2<T>& w) {
  // This is a general derivation.
  if (absLB(w.a - w.d) > 0) { // means axis(w) does meet (0,inf) orthogonally
    // normalized matrix so that w and sym(w) marg point is j
    SL2<T> normalized(w.a, sqrt(w.b * w.c), sqrt(w.b * w.c), w.d);
    if (g_debug && std::is_same<T, AJ>::value) {
      fprintf(stderr, "********** w_in_sym_search ***********\n");
      print_type("cosh_mu", cosh_move_j(normalized));
    }
    return absUB(cosh_move_j(normalized)) < g_cosh_sym_marg; 
  }
  return false; 
}

template<typename T>
inline const bool w_conj_in_sym_search(const SL2<T>& w,
    const Params<T>& p, const char g) {
  T t_sh_sq_hf_p = two_sinh_sqrd_half_perp_wg_zero_inf(w, p, g); 
  // does not intersect (0,inf) at a right angle
  // Note, if w(axis(g)) is (0, inf), then the margulis number
  // is at most len(g) since axis(g) and (0,inf) meet orthogonally
  if (absLB(t_sh_sq_hf_p + 1) > 0) {
    T ch_d = two_cosh_dist(t_sh_sq_hf_p) / 2;
    T coshl, cost;
    if (g == 'x') {
      coshl = p.coshlx;
      cost = p.costx;
    } else {
      coshl = p.coshly;
      cost = p.costy;
    }
    T cosh_mu = (ch_d * ch_d) * (coshl - cost) + cost;
    if (g_debug && std::is_same<T, AJ>::value) {
      fprintf(stderr, "********** w_conj_in_sym_search: %c ***********\n", g);
      print_type("two_sinh_sqrd_half_perp_wga_zero_inf", t_sh_sq_hf_p);
      print_type("cosh_dist_wga_zero_inf", ch_d);
      print_type("cosh_mu", cosh_mu);
    }
    if (absUB(cosh_mu) < g_cosh_sym_marg) {
      return true;
    }
  }
  return false;
}

template<typename T>
inline const bool w_conj_and_g_in_sym_search(const SL2<T>& w,
    const Params<T>& p, const char g) {
  T f_sh_sq_hf_p = four_sinh_sqrd_half_perp(w, p, g);
  if (absLB(f_sh_sq_hf_p) > 0 && absLB(f_sh_sq_hf_p + 4) > 0) {
    T ch_2r = four_cosh_dist(f_sh_sq_hf_p) / 4;
    // IMPORTANT: we assume that len(g) < g_sym_mu 
    // Thus, we only check that the distance to the
    // possible margulis point is small enough
    T tr;
    if (g == 'x') {
      tr = p.coshLx2;
    } else {
      tr = p.coshLy2;
    }
    T four_cosh_mu = ch_2r * (tr - 2) * (tr + 2) + abs_sqrd(tr); 
    if (g_debug && std::is_same<T, AJ>::value) {
      fprintf(stderr, "********** w_conj_in_sym_search: %c ***********\n", g);
      print_type("four_sinh_sqrd_half_perp_ga_wga", f_sh_sq_hf_p);
      print_type("cosh_dist_ga_wga", ch_2r);
      print_type("four_cosh_mu", four_cosh_mu);
    }
    if (absUB(four_cosh_mu) < g_cosh_sym_marg * 4) {
      return true;
    }
  }
  return false;
}

template<typename T>
inline const bool cant_fix_axis(const SL2<T>& w,
    const Params<T>& p, const char g) {
  T f_sh_sq_hf_p = four_sinh_sqrd_half_perp(w, p, g);
  if (g_debug && std::is_same<T, AJ>::value
      && absLB(f_sh_sq_hf_p) > 0 && absLB(f_sh_sq_hf_p + 4) > 0) {
    fprintf(stderr, "********** CANNOT FIX %c AXIS ***********\n", g);
    print_type(f_sh_sq_hf_p);
    fprintf(stderr, "LB values %f and %f\n", 
        absLB(f_sh_sq_hf_p), absLB(f_sh_sq_hf_p + 4));
    fprintf(stderr, "LB away from %d and %d\n",
        absLB(f_sh_sq_hf_p) > 0, absLB(f_sh_sq_hf_p + 4) > 0);
    fprintf(stderr, "*******************************\n");
  }
  return absLB(f_sh_sq_hf_p) > 0 && absLB(f_sh_sq_hf_p + 4) > 0; 
}

template<typename T>
inline const bool does_not_fix_zero_inf(const SL2<T>& w) {
  return (absLB(w.b) > 0 && absLB(w.d) > 0) || 
    (absLB(w.a) > 0 && absLB(w.c) > 0);
}

template<typename T>
inline const bool must_fix_axis(const SL2<T>& w,
    const Params<T>& p, const char g) {
  // The "must" part is only valid for AJ tests
  T ch_two_re_tube;
  if (g == 'x') {
    ch_two_re_tube = p.cosh2dx;
  } else {
    ch_two_re_tube = p.cosh2dy;
  }
  T diff = ch_two_re_tube * 4 - four_cosh_dist(w, p, g);
  if (g_debug && std::is_same<T, AJ>::value
      && strictly_pos(diff)) {
    fprintf(stderr, "********** MUST FIX %c AXIS ***********\n", g);
    print_SL2(w);
    print_type("4 cosh 2 re tube:", ch_two_re_tube * 4);
    print_type("4 cosh dist axis(g) and w(axis(g):",
        four_cosh_dist(w, p, g));
    print_type("diff:", diff);
    fprintf(stderr, "*******************************\n");
  }
  // We know that diff is away from zero and 
  // the diff should be conj symmetric, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
}

template<typename T>
inline const bool inside_var_nbd_g(const SL2<T>& w,
    const Params<T>& p, const char g) {
  // The second test may only work when g has trace close to +/- 2
  if (g_debug &&
      (absUB(jorgensen_wg(w, p, g)) < 1 ||
      absUB(jorgensen_gw(w, p, g)) < 1 ||
      must_fix_axis(w, p, g))) {
      fprintf(stderr, "UB Jwx %f, UB Jxw %f, must_fix %d\n",
          absUB(jorgensen_wg(w, p, g)),
          absUB(jorgensen_gw(w, p, g)),
          must_fix_axis(w, p, g));
  }
  return absUB(jorgensen_wg(w, p, g)) < 1 ||
      absUB(jorgensen_gw(w, p, g)) < 1 || must_fix_axis(w, p, g);
}

template<typename T>
inline const bool moves_y_axis_too_close_to_x(const SL2<T>& w,
    const Params<T>& p) {
  T diff = p.coshdxdy * 4 - four_cosh_dist_xwy(w, p);
  // We know that diff is away from zero and
  // the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  if (g_debug && std::is_same<T, AJ>::value && strictly_pos(diff)) {
    fprintf(stderr, "****************************************\n");
    fprintf(stderr, "MOVES Y TOO CLOSE TO X\n");
    print_SL2(w);
    print_type("4cosh(dx+dy):", p.coshdxdy * 4);
    print_type("4coshd(dist(x-axis, w(y-axis))):", 
        four_cosh_dist_xwy(w, p)); 
    print_type("diff:", diff);
    fprintf(stderr, "diff is positive: %d\n", strictly_pos(diff));
    fprintf(stderr, "****************************************\n");
  }
  return strictly_pos(diff);
}

template<typename T>
inline const bool moves_x_axis_too_close_to_y(const SL2<T>& w,
    const Params<T>& p) {
  return moves_y_axis_too_close_to_x(inverse(w), p);
}

template<typename T>
inline const bool moved_y_axis_not_x_axis(const SL2<T>& w, const Params<T>& p) {
  T f_sh_sq_hf_p = four_cosh_dist_xwy(w, p);
  return absLB(f_sh_sq_hf_p) > 0 && absLB(f_sh_sq_hf_p + 4) > 0; 
}

template<typename T>
inline const bool moved_x_axis_not_y_axis(const SL2<T>& w,
    const Params<T>& p) {
  return moved_y_axis_not_x_axis(inverse(w), p);
}

template<typename T>
inline bool margulis_smaller_than_xy(const SL2<T>& w1, const SL2<T>& w2, const Params<T>& p) {
  T diff = p.coshmu * 4 - four_cosh_margulis_simple(w1, w2).first;
  // We know that diff is away from zero and
  // the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
}

template<typename T>
inline bool move_less_than_marg(const SL2<T>& w, const Params<T>& p) {
  T diff = p.coshmu - cosh_move_j(w);
  return strictly_pos(diff); 
}

template<typename T>
inline bool non_cylic_power(const SL2<T>& w, const SL2<T>& g) {
  // Assume word fixes the same axis as x or y,
  // so it must live in a cyclic group with x or y.
  // Here we check that this is impossible in this box. Must use margulis
  // number to check cut off for roots of x or y
  SL2<T> commutator = g * w * inverse(w * g);
  if (g_debug && std::is_same<T, AJ>::value && not_identity(commutator)) {
    fprintf(stderr, "****************************************\n");
    fprintf(stderr, "NOT CYCLIC POWER\n");
    fprintf(stderr, "x or y\n");
    print_SL2(g);
    fprintf(stderr, "(x or y)^2\n");
    print_SL2(g * g);
    fprintf(stderr, "w\n");
    print_SL2(w);
    fprintf(stderr, "commutator\n");
    print_SL2(commutator);
    fprintf(stderr, "|b| == 0: %d, |c| == 0: %d, |a-1| == 0: %d, |d-1| == 0: %d, |a+1| == 0: %d, |d+1| == 0: %d\n", absLB(commutator.b) == 0, absLB(commutator.c) == 0, absLB(commutator.a-1) == 0,absLB(commutator.d-1) == 0, absLB(commutator.a+1) == 0, absLB(commutator.d+1) == 0);
    fprintf(stderr, "****************************************\n");
  }
  // TODO Test powers when coshmu > coshsdx + sinhsdx
  return not_identity(commutator); 
}

// Meyerhoff K Test
// We stop computing if we fail the test
#define MAX_MEYER 8
template<typename T>
bool meyerhoff_k_test(const T& ch_o, const T& cs_o,
    const T& four_cosh_tube_diam_UB) {
  // Assumed ch and cs are real valued jets for cosh(Re(L)) and cos(Im(L))
  T ch_prev = T(1);
  T cs_prev = T(1);
  T ch = ch_o;
  T cs = cs_o;
  T temp, four_cosh_tube_diam_LB;
  T meyer_k = T(1024); // arbitray large enough number
  double sqrt_of_2 = absLB(sqrt(Complex(2,0)));
  double sqrt_of_2_minus_one = absLB(sqrt(Complex(2,0)) - 1);
  int count = 0;
  while (absUB(ch) <= sqrt_of_2 && count < MAX_MEYER) {
    temp = ch - cs; 
    if (strictly_pos(meyer_k - temp) && absUB(temp) <= sqrt_of_2_minus_one) {
      meyer_k = temp;
      // See Meyerhoff paper on volume lowerbounds for hyperbolic 3-manifolds
      four_cosh_tube_diam_LB = sqrt(-(meyer_k * 32) + 16) / meyer_k;
      if (strictly_pos(four_cosh_tube_diam_LB - four_cosh_tube_diam_UB)) {
        if (g_debug) {
          fprintf(stderr,
              "Meyer k %f with 4 cosh tube diam LB %f and UB %f\n",
              absLB(meyer_k), absUB(four_cosh_tube_diam_LB),
              absLB(four_cosh_tube_diam_UB));
        }
        return true; // box can be killed
      }
    } 
    // Use Chebyshev recurrance relation
    temp = (ch_o * 2) * ch - ch_prev;
    ch_prev = ch;
    ch = temp;  
    temp = (cs_o * 2) * cs - cs_prev;
    cs_prev = cs;
    cs = temp;
    count +=1;
  }
  return false; // inconclusive
}

#define MAX_ROOTS 8
template<typename T>
T worst_primitive_cosh_re_len(const T& ch_o, const T& cs_o,
    const T& four_cosh_tube_diam_UB) {
  // Assumed ch and cs are real valued jets for cosh(Re(L)) and cos(Im(L))
  T ch_prev = ch_o;
  T cs_prev = cs_o;
  for (int i = 0; i < MAX_ROOTS; ++i) {
    T ch = sqrt((ch_prev + 1) / 2);
    // note, - pi <= Im(L) <= pi, so sign is +
    T cs = sqrt((cs_prev + 1) / 2);
    if (meyerhoff_k_test(ch, cs, four_cosh_tube_diam_UB)) {
      return ch_prev;
    }
    ch_prev = ch;
    cs_prev = cs;
  }
  // no luck
  T zero(0);
  return zero; 
}

template<typename T>
T cosh_marg_lower_bound(const T& two_sinh_r) {
  T s = two_sinh_r;
  T a8 = powT(s, 8) * (-0.0000014461700558); 
  T a7 = powT(s, 7) *   0.0000365880448817; 
  T a6 = powT(s, 6) * (-0.0003163830157272);
  T a5 = powT(s, 5) *   0.0005316504647188;
  T a4 = powT(s, 4) *   0.0086912125268823;
  T a3 = powT(s, 3) * (-0.061949675652791);
  T a2 = powT(s, 2) *   0.151649220047696;
  T a1 = s          * (-0.01513801009421);
  double a0 = 0.9999999;
  return ((a8 + (a1 + a0)) + (a4 + a5)) + ((a7 + a2) + (a6 + a3)); 
}

#define MAX_ID_SHIFT 5
template<typename T>
std::string proven_identity(std::string word, const Params<T>& p) {
  SL2<T> x = construct_x(p);
  SL2<T> y = construct_y(p);
  if (g_debug) {
    fprintf(stderr, "Testing proven identity for word: %s .\n", word.c_str());
  }
  SL2<T> w = construct_word(word, p);
  std::string new_word;
  if (y_power(word) > 0 && inside_var_nbd_g(w, p, 'x')) {
    T four_cosh_x_tube_UB = four_cosh_dist(y, p, 'x');
    T cosh_prim_re_len = worst_primitive_cosh_re_len(
        p.coshlx, p.costx, four_cosh_x_tube_UB); 
    for (auto s : {"x", "X"}) {
      new_word = x_strip(word);
      for (int i = 0; i < MAX_ID_SHIFT; ++i) {
        SL2<T> new_w = construct_word(new_word, p); // order matters
        T diff = cosh_prim_re_len * 4 - four_cosh_re_length(new_w);
        if (strictly_pos(diff)) {
          if (g_debug) {
            fprintf(stderr, "Found proven identity: %s .\n", new_word.c_str());
          }
          return new_word;
        }      
        new_word = s + new_word;
      }
    }
  }
  if (x_power(word) > 0 && inside_var_nbd_g(w, p, 'y')) {
    T four_cosh_y_tube_UB = four_cosh_dist(x, p, 'y');
    T cosh_prim_re_len = worst_primitive_cosh_re_len(
        p.coshly, p.costy, four_cosh_y_tube_UB); 
    for (auto s : {"y", "Y"}) {
      new_word = y_strip(word);
      for (int i = 0; i < MAX_ID_SHIFT; ++i) {
        SL2<T> new_w = construct_word(new_word, p); // order matters
        T diff = cosh_prim_re_len * 4 - four_cosh_re_length(new_w);
        if (strictly_pos(diff)) {
          if (g_debug) {
            fprintf(stderr, "Found proven identity: %s .\n", new_word.c_str());
          }
          return new_word;
        }      
        new_word = s + new_word;
      }
    }
  }
  return "";
}


#endif //_TestCollection_
