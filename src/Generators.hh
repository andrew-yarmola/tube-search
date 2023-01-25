#ifndef __Generators_h
#define __Generators_h
#include <math.h>
#include <unordered_map>
#include <string>
#include "SL2.hh"
#include "types.hh"

#define CACHE_WORD_LEN 9

template<typename T>
SL2<T> construct_x(const Params<T>& params) {
	return SL2<T>(params.coshLx2, params.expmdx * params.sinhLx2,
                params.expdx * params.sinhLx2, params.coshLx2);
};

template<typename T>
SL2<T> construct_y(const Params<T>& params) {
	return SL2<T>(params.coshLy2, params.expdyf * params.sinhLy2,
                params.expmdyf * params.sinhLy2, params.coshLy2);
};

template<typename T>
SL2<T> construct_word(std::string word, const Params<T>& params)
{
  std::unordered_map<std::string, SL2<T> > short_word_cache;
  return construct_word(word, params, short_word_cache);
/*  // TODO : add caching of powers if we are too slow
  // std::string recover = "";
  SL2<T> w; // identity
  SL2<T> x = construct_x(params);
  SL2<T> y = construct_y(params);

  char h;
  int x_pow = 0;
  int y_pow = 0;	
  std::string::reverse_iterator rit;
  for (rit = word.rbegin(); rit != word.rend(); ++rit) {
    h = *rit;
	  switch(h) {
      case 'x': ++x_pow; break;
      case 'X': --x_pow; break;
      case 'y': ++y_pow; break;
      case 'Y': --y_pow; break;
    }
    if (y_pow != 0 && x_pow != 0) {
      if (h == 'y' || h == 'Y' ) {
        w = pow(x, x_pow) * w;
        if (x_pow < 0) {
          recover = repeat("X",-x_pow) + recover;
        } else {
          recover = repeat("x",x_pow) + recover;
        }
        x_pow = 0;
      } else {
        w = pow(y, y_pow) * w;
        if (y_pow < 0) {
          recover = repeat("Y",-y_pow) + recover;
        } else {
          recover = repeat("y",y_pow) + recover;
        }
        y_pow = 0;
      }
    }
  }
  // Only one of these should be true
  if (x_pow != 0) { 
    w = pow(x, x_pow) * w;
    if (x_pow < 0) {
      recover = repeat("X",-x_pow) + recover;
    } else {
      recover = repeat("x",x_pow) + recover;
    }
  }
  if (y_pow != 0) { 
    w = pow(y, y_pow) * w;
    if (y_pow < 0) {
      recover = repeat("Y",-y_pow) + recover;
    } else {
      recover = repeat("y",y_pow) + recover;
    }
  }
  // fprintf(stderr, "Original %s vs %s\n", word.c_str(), recover.c_str());
	return w;
*/
};

template<typename T>
SL2<T> construct_word(std::string word, const Params<T>& params, std::unordered_map< std::string, SL2<T> >& word_cache) {
  if (word_cache.size() == 0) {
    SL2<T> x = construct_x(params);
    SL2<T> y = construct_y(params);
    word_cache.emplace("", SL2<T>());
    word_cache.emplace("x", x); 
    word_cache.emplace("y", y);
    word_cache.emplace("X", inverse(x)); 
    word_cache.emplace("Y", inverse(y));
  }
  // Check, just in case we've already seen this word
  auto lookup = word_cache.emplace(word, SL2<T>());
  if (lookup.second) { // not found
    word_cache.erase(lookup.first);
  } else {
    return lookup.first->second;
  }
  
  SL2<T> w_left = construct_word(word.substr(0, word.length()/2), params, word_cache);
  SL2<T> w_right = construct_word(word.substr(word.length()/2), params, word_cache);

  SL2<T> w = w_left * w_right;

  if (word.length() <= CACHE_WORD_LEN) {
    word_cache.emplace(word, w);
  }

  return w;
}


#endif // __Generators_h
