#include <list>
#include <map>
#include <queue>
#include "TubeSearch.hh"
#include "TestCollection.hh"
#include "CanonicalName.hh"
#include "IsomH3.hh"
#include "Params.hh"
#include "roundoff.h"

using namespace std;

namespace TubeSearchImpl {

  inline bool pos_less(double a, double b) {
    return a > 0 && a < b;
  }

  inline bool pos_neq(double a, double b) {
    return a > 0 && b > 0 && a != b;
  }

  struct Word {
    string name;
    SL2<Complex> matrix;
    void compute_vals(const Params<Complex>& params);
    double four_cosh_x_x = -1;
    double four_cosh_x_y = -1;
    double four_cosh_y_y = -1;
    double four_cosh_y_x = -1;
    double four_cosh_re_len = -1;
    double cosh_move = -1;
    double jorg_x = -1;
    double jorg_y = -1;
    bool operator<(const Word& b) const {
      if (x_power(name) > 0 && y_power(name) > 0 && x_power(b.name) > 0 && y_power(b.name) > 0) {
        if (pos_neq(four_cosh_x_x, b.four_cosh_x_x)) {
          return four_cosh_x_x > b.four_cosh_x_x;
        }
        if (pos_neq(four_cosh_y_y, b.four_cosh_y_y)) {
          return four_cosh_y_y > b.four_cosh_y_y;
        }
        if (pos_neq(four_cosh_x_y, b.four_cosh_x_y)) {
          return four_cosh_x_y > b.four_cosh_x_y;
        }
        if (pos_neq(four_cosh_y_x, b.four_cosh_y_x)) {
          return four_cosh_y_x > b.four_cosh_y_x;
        }
        if (pos_neq(cosh_move, b.cosh_move)) {
          return cosh_move > b.cosh_move;
        }
        if (pos_neq(jorg_x, b.jorg_x)) {
          return jorg_x > b.jorg_x;
        }
        if (pos_neq(jorg_y, b.jorg_y)) {
          return jorg_y > b.jorg_y;
        }
        if (pos_neq(four_cosh_re_len, b.four_cosh_re_len)) {
          return four_cosh_re_len > b.four_cosh_re_len;
        }
      } else {
        if (x_power(name) > 0 && y_power(name) > 0) {
          return false;
        }
        if (x_power(b.name) > 0 && y_power(b.name) > 0) {
          return true;
        }
      }
      return name > b.name;
    }
  };

  struct WordPair {
    WordPair(Word& first, Word& second, const Params<Complex>& params);
    Word first;
    Word second;
    double cosh_two_ortho = -1;
    double four_cosh_marg = -1;
    double jorgensen_words = -1;
  };

  struct PairComp : public binary_function<const WordPair&, const WordPair&, bool>
  {
    bool operator() (const WordPair& a, const WordPair& b)
    {
      if (pos_neq(a.cosh_two_ortho, b.cosh_two_ortho)) {
        return a.cosh_two_ortho < b.cosh_two_ortho;
      }
      if (pos_neq(a.four_cosh_marg, b.four_cosh_marg)) {
        return a.four_cosh_marg < b.four_cosh_marg;
      }
      if (pos_neq(a.jorgensen_words, b.jorgensen_words)) {
        return a.jorgensen_words < b.jorgensen_words;
      }
      if (a.first.name != b.first.name) {
        return a.first < b.first;
      }
      return a.second < b.second;
    }
  };

  string name_inverse(string name) {
    string inv = name;
    reverse(inv.begin(), inv.end());
    string::size_type pos;
    for (pos = 0; pos < inv.size(); ++pos) {
      char c = inv[pos];
      if (c >= 'a' && c <= 'z')
        inv[pos] += 'A' - 'a';
      else
        inv[pos] -= 'A' - 'a';
    }
    return inv;
  }

  Word inverse(Word& word, const Params<Complex> params)
  {
    Word inv;
    inv.name = name_inverse(word.name);
    inv.matrix = inverse(word.matrix);
    inv.compute_vals(params);
    return inv;
  }

  struct TubeSearch {
    TubeSearch(Params<Complex> params_) :params(params_), m_found_good_pair(false) {
      x_word.name = "x";
      x_word.matrix = construct_x(params);
      x_word.compute_vals(params);
      y_word.name = "y";
      y_word.matrix = construct_y(params);
      y_word.compute_vals(params);
      X_word = inverse(x_word, params);
      Y_word = inverse(y_word, params);
    }

    ~TubeSearch() {}
    bool push_word(Word word);
    bool push_word(string name);
    WordPair find_pair();
    bool found_good_pair() { return m_found_good_pair; }
    void add_relator(string word);

    const Params<Complex> params;
    CanonicalName canonical_name;
    Word x_word;
    Word y_word;
    Word X_word;
    Word Y_word;
    bool m_found_good_pair;
    vector<Word> word_lookup;
    priority_queue<Word> words;
    priority_queue<WordPair, vector<WordPair>, PairComp> pairs;
  };

  Word product(const Word& lhs, const Word& rhs, const Params<Complex>& params, CanonicalName& canonical_name)
  {
    Word result;
    result.name = canonical_name.get_canonical_name(lhs.name + rhs.name);
    if (result.name == lhs.name + rhs.name) {
      result.matrix = lhs.matrix * rhs.matrix;
    } else {
      result.matrix = construct_word(result.name, params);
    }
    result.compute_vals(params);
    return result;
  }

  void Word::compute_vals(const Params<Complex>& params)
  {
    four_cosh_re_len = absUB(four_cosh_re_length(matrix));
    if (y_power(name) > 0) {
      four_cosh_x_x = absUB(four_cosh_dist(matrix, params, 'x')); 
      four_cosh_y_x = absUB(four_cosh_dist_xwy(matrix, params)); 
      jorg_x = min(absUB(jorgensen_gw(matrix, params, 'x')),
          absUB(jorgensen_wg(matrix, params, 'x'))); 
      if (x_power(name) > 0) {
        cosh_move = absUB(cosh_move_j(matrix));
      }
    }
    if (x_power(name) > 0) {
      four_cosh_y_y = absUB(four_cosh_dist(matrix, params, 'y')); 
      four_cosh_x_y = absUB(four_cosh_dist_xwy(inverse(matrix), params)); 
      jorg_y = min(absUB(jorgensen_gw(matrix, params, 'y')),
          absUB(jorgensen_wg(matrix, params, 'y'))); 
    }
  }

  WordPair::WordPair(Word& first, Word& second, const Params<Complex>& params):first(first), second(second)
  {
    //  fprintf(stderr, "Building word pair\n");
    //  fprintf(stderr, "a_1 = %f + %f i, c_1 = %f + %f i\n a_2 =  %f + %f i, c_2 = %f + %f i\n",
    //                   first.matrix.a.re, first.matrix.a.im, first.matrix.c.re, first.matrix.c.im,
    //                   second.matrix.a.re, second.matrix.a.im, second.matrix.c.re, second.matrix.c.im); 
    if (first.name.length() != 0 && second.name.length() != 0) {
      cosh_two_ortho = absUB(cosh_2_re_perp(first.matrix, second.matrix));  
      jorgensen_words = absUB(jorgensen(first.matrix, second.matrix));
      four_cosh_marg = absUB(four_cosh_margulis_simple(first.matrix, second.matrix).first);
    }
    //    fprintf(stderr, "Built pair (%s,%s)\n",
    //      first.name.c_str(), second.name.c_str());
  }

  bool TubeSearch::push_word(Word word)
  {
    if (word.name == "") { return false; }
    for (auto it = word_lookup.begin(); it != word_lookup.end(); ++it) {
      if (it->name == word.name) {
        return false;
      }
    }
    // We assume compute_vals has run on the word
    // fprintf(stderr, "Word: %s\n", word.name.c_str());
    bool pushed = false;

    if (pos_less(word.four_cosh_x_x, absUB(params.coshdxdy * 6.5))   ||  
        pos_less(word.four_cosh_y_y, absUB(params.coshdxdy * 6.5))   ||
        pos_less(word.four_cosh_x_y, absUB(params.coshdxdy * 6.5))  ||
        pos_less(word.four_cosh_y_x, absUB(params.coshdxdy * 6.5))  ||
        pos_less(word.cosh_move, absUB(params.coshmu * 2.5))    ||
        pos_less(word.jorg_x, 3)                                ||
        pos_less(word.jorg_y, 3)                                ||
        pos_less(word.four_cosh_re_len, absUB(params.coshlx * 6.5)))
    {
      word_lookup.push_back(word);
      words.push(word);
      pushed = true;
      //fprintf(stderr, "Pushed word: %s\n", word.name.c_str());
    }

    if (x_power(word.name) == 0 || y_power(word.name) == 0) {
      return pushed;
    }

    //    return pushed;

    for (auto it = word_lookup.begin(); it != word_lookup.end(); ++it) {
      //    fprintf(stderr,"Trying to print\n");
      //    fprintf(stderr, "First word %s, Second word %s\n", it->name.c_str(), word.name.c_str());
      if (word.name == it->name || word.name == name_inverse(it->name) || x_power(it->name) == 0 || y_power(it->name) == 0) {
        continue;
      }
      WordPair pair(word, *it, params);
      if (pos_less(pair.jorgensen_words, 1.5)                        ||
          pos_less(pair.four_cosh_marg, absUB(params.coshmu * 4.1))    ||
          pos_less(pair.cosh_two_ortho, absUB(params.coshdxdy * 1.2)))
      {
        pairs.push(pair);
        // pushed = true;
        // fprintf(stderr,"Pushed pair (%s,%s)\n", pair.first.name.c_str(), pair.second.name.c_str());
      }
    }
    return pushed;
    //  fprintf(stderr,"Done with push_word\n");
  }

  bool TubeSearch::push_word(string name)
  {
    name = canonical_name.get_canonical_name(name);
    if (name == "") return false;
    Word w;
    w.name = name;
    w.matrix = construct_word(w.name, params);
    w.compute_vals(params);
    return push_word(w);
    //	Word wInv = inverse(w);
    //	find_names(wInv);
    //	push_word(wInv);
    //  fprintf(stderr, "pushed words %s(%s) and %s(%s)\n",
    //  w.name.c_str(), w.nameClass.c_str(),
    //  wInv.name.c_str(), wInv.nameClass.c_str());
  }

#define SEARCH_DEPTH 30
#define POW_DEPTH 6

  WordPair TubeSearch::find_pair()
  {
    int loop_count = 0;
    while ((words.size() > 0 || pairs.size() > 0) && loop_count < SEARCH_DEPTH) {
      if (words.size() > 0) {
        Word w(words.top());
        // fprintf(stderr, "Smallest word %s\n", w.name.c_str());
        //        if (w.name == "Xy") {
        //          fprintf(stderr, "    4cosh(ax, wax): %f < %f\n    4cosh(ay, way): %f < %f\n    4cosh(ax, way): %f < %f\n    4cosh(ay, wax): %f < %f\n    cosh_move: %f < %f\n    jorg_x: %f < 1.0\n    jorg_y: %f < 1.0\n", w.four_cosh_x_x, absLB(params.coshdxdy * 4), w.four_cosh_y_y, absLB(params.coshdxdy * 4), w.four_cosh_x_y, absLB(params.coshdxdy * 4), w.four_cosh_y_x, absLB(params.coshdxdy * 4), w.cosh_move, absLB(params.coshmu), w.jorg_x, w.jorg_y);
        //        } 
        if (y_power(w.name) > 0 && x_power(w.name) > 0 &&
            (pos_less(w.four_cosh_x_x, absLB(params.coshdxdy * 3.999999))   ||  
             pos_less(w.four_cosh_y_y, absLB(params.coshdxdy * 3.999999))   ||
             pos_less(w.four_cosh_x_y, absLB(params.coshdxdy * 3.999999))  ||
             pos_less(w.four_cosh_y_x, absLB(params.coshdxdy * 3.999999))  ||
             pos_less(w.cosh_move, absLB(params.coshmu * 0.999999))        ||
             pos_less(w.jorg_x, 0.9999999)                                 ||
             pos_less(w.jorg_y, 0.9999999)                                ||
             pos_less(w.four_cosh_re_len, absLB(params.coshlx * 3.9999999))))
        {
          m_found_good_pair = true;
          Word none;
          none.name = "";
          WordPair result(w, none, params);
          // fprintf(stderr, "Found word %s\n", w.name.c_str());
          // fprintf(stderr, "    4cosh(ax, wax): %f < %f\n    4cosh(ay, way): %f < %f\n    4cosh(ax, way): %f < %f\n    4cosh(ay, wax): %f < %f\n    cosh_move: %f < %f\n    jorg_x: %f < 1.0\n    jorg_y: %f < 1.0\n", w.four_cosh_x_x, absLB(params.coshdxdy * 4), w.four_cosh_y_y, absLB(params.coshdxdy * 4), w.four_cosh_x_y, absLB(params.coshdxdy * 4), w.four_cosh_y_x, absLB(params.coshdxdy * 4), w.cosh_move, absLB(params.coshmu), w.jorg_x, w.jorg_y);
          return result;
        }
        words.pop();
        if (w.name[0] == 'y' || w.name[0] == 'Y') {
          int x_pow = 0;
          Word x_mult = product(x_word, w, params, canonical_name);
          while(push_word(x_mult) && x_pow < POW_DEPTH) {
            x_mult = product(x_word, x_mult, params, canonical_name);
            x_pow++;
          }
          int X_pow = 0;
          Word X_mult = product(X_word, w, params, canonical_name);
          while(push_word(X_mult) && X_pow < POW_DEPTH) {
            X_mult = product(X_word, X_mult, params, canonical_name);
            X_pow++;
          }
        } else {
          int y_pow = 0;
          Word y_mult = product(y_word, w, params, canonical_name);
          while(push_word(y_mult) && y_pow < POW_DEPTH) {
            y_mult = product(y_word, y_mult, params, canonical_name);
            y_pow++;
          }
          int Y_pow = 0;
          Word Y_mult = product(Y_word, w, params, canonical_name);
          while(push_word(Y_mult) && Y_pow < POW_DEPTH) {
            Y_mult = product(Y_word, Y_mult, params, canonical_name);
            Y_pow++;
          }
        }
      }
      if (pairs.size() > 0) { 
        WordPair p(pairs.top());
        // fprintf(stderr, "Smallest pair (%s, %s)\n", p.first.name.c_str(), p.second.name.c_str());
        if (pos_less(p.jorgensen_words, 0.999999)                ||
            pos_less(p.four_cosh_marg, absLB(params.coshmu * 3.999999))) 
          //pos_less(p.cosh_two_ortho, absLB(params.coshdxdy)))
        {
          m_found_good_pair = true;
          pairs.pop();
          // fprintf(stderr, "Found pair (%s, %s)\n", p.first.name.c_str(), p.second.name.c_str());
          // fprintf(stderr, "    jorgensen_words: %f < 1.0\n    four_cosh_marg: %f < %f\n    cosh_two_ortho: %f < %f\n", p.jorgensen_words, p.four_cosh_marg, absLB(params.coshmu * 4), p.cosh_two_ortho, absLB(params.coshdxdy));
          return p;
        }
        pairs.pop();
        push_word(product(inverse(p.first, params), p.second, params, canonical_name));
        //push_word(product(inverse(p.second, params), p.first, params, canonical_name));
        //push_word(product(p.first, inverse(p.second, params), params, canonical_name));
        //push_word(product(p.first, p.second, params, canonical_name));
      }
      loop_count++;
    }
    Word w;
    w.name = "";
    return WordPair(w, w, params);
  }

  void TubeSearch::add_relator(string w)
  {
    canonical_name.add_relator(w);
  }
}

#define MAX_ATTEMPTS 2

vector< word_pair > find_pairs(Params<Complex> center, vector<string> seed_words,
    int num_words, int max_length, vector<string> relators)
{
  TubeSearchImpl::TubeSearch search(center);
  search.add_relator("xX"); 
  search.add_relator("Xx"); 
  search.add_relator("yY"); 
  search.add_relator("Yy"); 
  for (int i = 0; i < relators.size(); ++i)
    if (relators[i].size() > 0)
      search.add_relator(relators[i]);
  search.push_word("x");
  search.push_word("y");
  search.push_word("X");
  search.push_word("Y");
  for (int i = 0; i < seed_words.size(); ++i)
    if (seed_words[i].size() > 0)
      search.push_word(seed_words[i]);
  set< word_pair > found_pairs;
  int attempts = 0;
  while (attempts < MAX_ATTEMPTS && (num_words > int(found_pairs.size()) || (-num_words > int(found_pairs.size()) && !search.found_good_pair()))) {
    //fprintf(stderr, "While loop: words size %d\n", (int) found_pairs.size());
    TubeSearchImpl::WordPair p = search.find_pair();
    //if (search.found_good_pair()) {
    //   printf("FOUND PAIR\n");
    //}
    attempts++;
    if (p.first.name.length() == 0 || (p.first.name.length() > max_length && p.second.name.length() > max_length)) {
      search.m_found_good_pair = false;
      continue;
    }
    word_pair P(p.first.name, p.second.name);
    found_pairs.insert(P);
  }
  vector< word_pair > v(found_pairs.begin(), found_pairs.end());
  return v;
}

typedef struct {
  Complex p; // attracting
  Complex m; // repelling
  SL2<Complex> gamma;
  string word;
} axis; 

//inline double sqnorm(const Complex &x) {
//    return x.re * x.re + x.im * x.im;
//}

//inline const Complex conj(const Complex &x) {
//    return Complex(x.re, -x.im);
//}

double cosh_re_orth_LB(const axis &u, const axis &v)
{
  Complex n1 = (u.p - v.m) * (u.m - v.p); 
  Complex n2 = (u.m - v.m) * (u.p - v.p);
  Complex d  = (u.m - u.p) * (v.m - v.p);
  return (absLB(n1) + absLB(n2))/absUB(d);    
}

double g_cut_ratio = 3.0;
double g_eps = pow(2,-100);

string word_inverse(const string &word)
{
  string inv = word;
  reverse(inv.begin(), inv.end());
  for_each(inv.begin(), inv.end(), [](char & c){
      c = islower(c) ? toupper(c) : tolower(c);
      });
  return inv;
}

inline void move(axis& a, const string& word, const SL2<Complex>& gamma) {
  a.p = mobius(gamma, a.p);
  a.m = mobius(gamma, a.m);
  a.word = word + a.word;
  a.gamma = gamma * a.gamma;
}  

#define MAX_SEEN_AGAIN 128
#define MAX_SHIFT 7
#define MAX_TOTAL 1000000

vector<string> find_words_tubes(const axis &to_move, bool x_is_shifter,
    bool x_is_mover, const Params<Complex> &params,
    double cosh_ortho_bound, int num_words, int max_levels,
    const vector<string>& relators, const map<string, int>& seen)
{
  vector<string> new_words;
  vector<axis> level_zero;
  level_zero.push_back(to_move);
  map< int, vector<axis> > axes;
  axes[0] = level_zero;
  // Generate new axes
  int d = 0;
  int seen_count = 0;
  int total = 0;
  SL2<Complex> x = construct_x(params);
  SL2<Complex> X = inverse(x);
  SL2<Complex> y = construct_y(params);
  SL2<Complex> Y = inverse(y);
  SL2<Complex> shifter;
  SL2<Complex> shifter_inv;
  axis fixed;
  string shift_word;
  string shift_word_inv;
  if (x_is_shifter) {
    fixed.p =  params.expmdx;
    fixed.m = -params.expmdx;
    shifter = x;
    shifter_inv = X;
    shift_word = "x";
    shift_word_inv = "X";
  } else {
    fixed.p =  params.expdyf;
    fixed.m = -params.expdyf;
    shifter = y;
    shifter_inv = Y;
    shift_word = "y";
    shift_word_inv = "Y";
  }
  map< string, SL2<Complex> > valid;
  while (d < max_levels) {
    vector<axis> level;
    axes[d+1] = level;
    for (const auto &a : axes[d]) {
      string first = a.word.substr(0,1);
      if (x_is_mover) {
        valid = {{"x", x}, {"X", X}};
      } else {
        valid = {{"y", y}, {"Y", Y}};
      }
      valid.erase(word_inverse(first));
      // Apply X,x or Y,y if possible
      for (const auto &h : valid) {
        axis h_moved = a;
        move(h_moved, h.first, h.second);
        double c_re_orth = cosh_re_orth_LB(fixed, h_moved);
        // fprintf(stderr, "Axis with word %s has %f distance vs %f\n",
        // h_moved.word.c_str(), c_re_orth, cosh_ortho_bound);
        double c_move_j = absUB(cosh_move_j(h_moved.gamma));
        Complex trace = h_moved.gamma.a + h_moved.gamma.d;
        double f_cosh_re_len = absUB(four_cosh_re_length(h_moved.gamma));

        if (c_move_j < absLB(0.98 * params.coshmu) ||
            (absUB(trace - 2) < 0.1 || absUB(trace + 2) < 0.1) || 
            c_re_orth < 0.98 * cosh_ortho_bound ||
            f_cosh_re_len < absLB(params.coshlx * 3.99)) {
          if (find(relators.begin(), relators.end(), h_moved.word) == relators.end()) {
            if (seen.find(h_moved.word) == seen.end()) {
              /*if (c_re_orth < 0.98 * cosh_ortho_bound) {
                fprintf(stderr, "Axis with word %s has %f distance vs %f\n",
                h_moved.word.c_str(), c_re_orth, cosh_ortho_bound);
                } else if (c_move_j < absLB(0.98 * params.coshmu)) {
                fprintf(stderr, "Word %s has move j %f vs %f\n",
                h_moved.word.c_str(), c_move_j, absLB(params.coshmu));
                }*/
              new_words.push_back(h_moved.word);
            } else {
              seen_count += 1;
            }
            if (new_words.size() >= num_words || seen_count >= MAX_SEEN_AGAIN) {
              return new_words;
            }
          }
        } else if (c_re_orth > 2.5 * cosh_ortho_bound) {
          continue;
        } else  {
          axes[d+1].push_back(h_moved);
          total += 1;
          axis shifted = h_moved;
          axis shifted_inv = h_moved;
          int shift_count = 0;
          while  (shift_count < MAX_SHIFT) {
            move(shifted, shift_word, shifter); 
            move(shifted_inv, shift_word_inv, shifter_inv); 
            axes[d+1].push_back(shifted);
            axes[d+1].push_back(shifted_inv);
            shift_count += 1 ;
            total += 2;
          } 
        }
        if (total > MAX_TOTAL) {
          return new_words;
        }
      }
    }
    d += 1;
  }
  return new_words;
}


vector< word_pair > find_words_v2(const Params<Complex>& params, int num_words, int max_move_len,
    const vector<string>& relators, const map<string, int>& seen)
{
  axis x_axis, y_axis;
  x_axis.p = params.expmdx;
  x_axis.m = -params.expmdx;
  y_axis.p = params.expdyf;
  y_axis.m = -params.expdyf;
  vector<string> move_x_words = find_words_tubes(x_axis, true, false, params,
      absUB(params.coshdxdy), num_words, max_move_len, relators, seen);
  vector<string> move_y_words = find_words_tubes(y_axis, false, true, params,
      absUB(params.coshdxdy), num_words, max_move_len, relators, seen);
  vector<string> move_xy_words = find_words_tubes(y_axis, true, false, params,
      absUB(params.coshdxdy), num_words, max_move_len, relators, seen);
  vector<string> move_yx_words = find_words_tubes(x_axis, false, true, params,
      absUB(params.coshdxdy), num_words, max_move_len, relators, seen);

  set<string> new_words;
  new_words.insert(move_x_words.begin(), move_x_words.end());
  new_words.insert(move_y_words.begin(), move_y_words.end());
  new_words.insert(move_xy_words.begin(), move_xy_words.end());
  new_words.insert(move_yx_words.begin(), move_yx_words.end());

  vector<word_pair> result;
  for (const auto &w : new_words) {
    word_pair wp(w, "");
    result.push_back(wp);
  }    
  return result;
}
