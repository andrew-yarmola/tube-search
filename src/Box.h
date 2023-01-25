#ifndef __Box_h
#define __Box_h
#include "types.hh"
#include "SL2.hh"
#include "AJ.h"
#include "QuasiRelators.h"
#include <unordered_map>

#define DIM 6
#define SCL 8   
// initial box radii are 2^(18/6), 2^(17/6), ..., 2^(13/6). The last is > 4.49

struct Box {
  Box();
  std::string name;
  std::string desc();
  QuasiRelators qr;
  Box child(int dir) const;
  Params<Complex> center() const { return _center; }
  // returns all values closer to 0 than in box or 0 if box overlaps
	// TODO Params<Complex> nearer() const { return _nearer; }
  Params<AJ> cover() const { return _cover; }
  SL2<Complex> x_center() const { return _x_center; }
  SL2<Complex> y_center() const { return _y_center; }
  SL2<AJ> x_cover() const { return _x_cover; }
  SL2<AJ> y_cover() const { return _y_cover; }
  std::unordered_map<std::string, SL2<AJ> > short_words_cache;
private:
  int pos;
  double center_digits[DIM];
  double size_digits[DIM];
  double box_center[DIM];
  double box_size[DIM];
  Params<Complex> _center;
  // Params<Complex> _nearer;
  Params<AJ> _cover;
  void compute_center_and_size();
  void compute_cover();
  // void compute_nearer();
  SL2<Complex> _x_center;
  SL2<Complex> _y_center;
  SL2<AJ> _x_cover;
  SL2<AJ> _y_cover;
};

Box get_box(std::string code);

#endif // __Box_h
