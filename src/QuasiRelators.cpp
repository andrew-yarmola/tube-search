#include <set>
#include <vector>
#include <algorithm>
#include "QuasiRelators.h"
#include "types.hh"
using namespace std;

string QuasiRelators::get_name(string w)
{
	NameStore::iterator it = names.find(w);
	if (it == names.end()) {
		add_quasi_relator(w);
		it = names.find(w);
	}
	return it->second;
}

void QuasiRelators::add_quasi_relator(string w)
{
	if (names.find(w) != names.end())
		return;
  // TODO: see if any kind of uniqueness is possible (probably not, since mult not commutative)
  names[w] = w;
	name_vector.push_back(w);
}

string QuasiRelators::min_pow_desc()
{
  vector<string> qrs(name_vector);
  sort(qrs.begin(), qrs.end(), x_power_sort);
  int min_power = x_power(qrs.front());
	string buf;
	for (vector<string>::iterator it = qrs.begin(); it != qrs.end(); ++it) {
    if (x_power(*it) > min_power) {
      break;
    }
		if (!buf.empty()) {
      buf += ",";
    } 
		buf += *it;
  }
	return buf;
}

string QuasiRelators::desc()
{
	string buf;
	for (vector<string>::iterator it = name_vector.begin(); it != name_vector.end(); ++it) {
		if (!buf.empty())
			buf += ",";
		buf += *it;
	}
	return buf;
}

std::vector<std::string> all_words();
std::vector<std::string> word_classes();

vector<string> QuasiRelators::word_classes()
{
	return name_vector;
}

vector<string> QuasiRelators::all_words()
{
	vector<string> result;
	for (NameStore::iterator it = names.begin(); it != names.end(); ++it) {
		result.push_back(it->first);
	}
	return result;
}

bool QuasiRelators::is_quasi_relator(std::string w)
{
	return names.find(w) != names.end();
}

string QuasiRelators::inverse(string w)
{
  reverse(w.begin(), w.end());
  string::size_type pos;
  for (pos = 0; pos < w.size(); ++pos) {
    char c = w[pos];
    if (c >= 'a' && c <= 'z')
      w[pos] += 'A' - 'a';
    else
      w[pos] -= 'A' - 'a';
  }
  return w;
}
