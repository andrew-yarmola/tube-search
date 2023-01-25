#include "types.hh"
#include "roundoff.h"

using namespace std;

double pi = atan(1.0) * 4;
double twopi = atan(1.0) * 8;

double absUB(const Complex& x) {
  return (1+2*EPS)*hypot(x.real(), x.imag());
}

double absLB(const Complex& x) {
  return (1-2*EPS)*hypot(x.real(), x.imag());
}

void split_string(const string &str, const string &delims, vector<string> &out) {
  size_t start;
  size_t end = 0;
  while ((start = str.find_first_not_of(delims, end)) != string::npos) {
    end = str.find_first_of(delims, start);
    out.push_back(str.substr(start, end - start));
  }
}

Complex shift_imag_around_zero(const Complex &a) {
  double im = a.imag();
  while (im < -pi) { im += twopi; }
  while (im >  pi) { im -= twopi; }
  return Complex(a.real(), im);
}

Complex shift_imag_around_pi(const Complex &a) {
  double im = a.imag();
  while (im <     0) { im += twopi; }
  while (im > twopi) { im -= twopi; }
  return Complex(a.real(), im);
}

Complex parse_complex(const string &complex_str) {
  vector<string> parts;
  // printf("complex: %s\n", complex_str.c_str());
  split_string(complex_str, "(-+j)", parts);
  if (complex_str.find("-") != string::npos) {
    return Complex(stod(parts[0]), -stod(parts[1]));
  } else {
    return Complex(stod(parts[0]), stod(parts[1]));
  }
}

string repeat(string s, int n) {
  string t = "";
  if (n <= 0) { return t; }
  string r = s;
  while (n > 1) {
    if (n & 1) { // n odd
      t += r;
    }
    r += r;
    n /= 2; // int division
  } 
  return t+r;
}

string cyclic_strip(string w) {
  size_t first = -1;
  size_t last = w.length();
  for (size_t c = 1; c < w.length(); ++c) {
    if ((w[first + c] == 'x' &&  w[last - c] == 'X') ||
        (w[first + c] == 'y' &&  w[last - c] == 'Y') ||
        (w[first + c] == 'X' &&  w[last - c] == 'x') ||
        (w[first + c] == 'Y' &&  w[last - c] == 'y')) {
      continue; 
    } else {
      first = first + c - 1; 
      last = first - c + 1; 
      break;
    }
  }
  return w.substr(first + 1, last - first - 1);
}

string x_strip(string w) {
  size_t first = -1;
  for (size_t p = 0; p < w.length(); ++p) {
    if (w[p] == 'x' || w[p] == 'X') {
      first = p; 
    } else {
      break;
    }
  }
  size_t last = w.length(); 
  for (size_t p = w.length() - 1; p >= 0; --p) {
    if (w[p] == 'x' || w[p] == 'X') {
      last = p; 
    } else {
      break;
    }
  }
  return w.substr(first+1, last - first - 1);
}

string y_strip(string w) {
  size_t first = -1;
  for (size_t p = 0; p < w.length(); ++p) {
    if (w[p] == 'y' || w[p] == 'Y') {
      first = p; 
    } else {
      break;
    }
  }
  size_t last = w.length();
  for (size_t p = w.length() - 1; p >= 0; --p) {
    if (w[p] == 'y' || w[p] == 'Y') {
      last = p; 
    } else {
      break;
    }
  }
  return w.substr(first+1, last - first - 1);
}

string x_rstrip(string w) {
  size_t last = w.length(); 
  for (size_t p = w.length() - 1; p >= 0; --p) {
    if (w[p] == 'x' || w[p] == 'X') {
      last = p; 
    } else {
      break;
    }
  }
  return w.substr(0, last);
}

string y_rstrip(string w) {
  size_t last = w.length();
  for (size_t p = w.length() - 1; p >= 0; --p) {
    if (w[p] == 'y' || w[p] == 'Y') {
      last = p; 
    } else {
      break;
    }
  }
  return w.substr(0, last);
}


int syllables(string w) {
  int count = 0;
  char cur = 'z'; // any char not in list
  for (string::size_type p = 0; p < w.size(); ++p) {
      if (w[p] != tolower(cur) && w[p] != toupper(cur)) {
        ++count;
        cur = w[p];
      }
  }
  return count;
} 

int x_power(string w) {
  int count = 0;
  for (string::size_type p = 0; p < w.size(); ++p) {
      if (w[p] == 'x' || w[p] == 'X') ++count;
  }
  return count;
} 

int y_power(string w) {
  int count = 0;
  for (string::size_type p = 0; p < w.size(); ++p) {
      if (w[p] == 'y' || w[p] == 'Y') ++count;
  }
  return count;
} 

bool x_power_sort(string a, string b) { return x_power(a) < x_power(b); }

bool y_power_sort(string a, string b) { return y_power(a) < y_power(b); }

template<>
void print_type<const Complex>(const Complex& x) {
  fprintf(stderr, "%f + %f I\n", x.real(), x.imag());
}

template<>
void print_type<Complex>(Complex& x) {
  print_type((const Complex) x);
}

template<>
void print_center<Complex>(const Complex& x) {
  print_type((const Complex) x);
}

template<>
bool sort_comp<Complex>(const Complex& a, const Complex& b) {
  return absUB(a) < absUB(b);
}

string double_to_hex(double x)
{
  char buf[100];
  union {
    double d;
    long l[2];
  } u;
  u.d = x;
  sprintf(buf, "'%0lx', '%0lx' (%.18f)", u.l[0], u.l[1], u.d);
  return string(buf);
}
