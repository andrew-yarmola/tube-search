#include "CanonicalName.hh"
#include <set>
#include <algorithm>
#include <stdio.h>

using namespace std;

map<string, CanonicalName*> CanonicalName::cache;

int x_pow(string w) {
  int count = 0;
  for (string::size_type p = 0; p < w.size(); ++p) {
      if (w[p] == 'x' || w[p] == 'X') ++count;
  }
  return count;
}; 

int y_pow(string w) {
  int count = 0;
  for (std::string::size_type p = 0; p < w.size(); ++p) {
      if (w[p] == 'y' || w[p] == 'Y') ++count;
  }
  return count;
}; 


char invertLetter(char c)
{
	switch (c) {
		case 'x': return 'X';
		case 'X': return 'x';
		case 'y': return 'Y';
		case 'Y': return 'y';
		default: return '?';
	}
}

string CanonicalName::inverse(string s)
{
	reverse(s.begin(), s.end());
	for (string::size_type sp = 0; sp < s.size(); ++sp)
		s[sp] = invertLetter(s[sp]);
	return s;
}

CanonicalName::CanonicalName()
{
	impl = 0;
}

void CanonicalName::initImpl()
{
	string implName;
	for (vector<string>::iterator it = relators.begin(); it != relators.end(); ++it) {
		implName += *it + ",";
	}
	impl = cache[implName];
	if (!impl) {
		impl = new CanonicalName();
		cache[implName] = impl;
		impl->initSubstitutions();
		for (vector<string>::iterator it = relators.begin(); it != relators.end(); ++it) {
			impl->add_relator_internal(*it);
		}
	}
}


void CanonicalName::initSubstitutions()
{
	substitutions.push_back(Substitution("xX", ""));
	substitutions.push_back(Substitution("Xx", ""));
	substitutions.push_back(Substitution("yY", ""));
	substitutions.push_back(Substitution("Yy", ""));
}

void CanonicalName::add_substitution(string& a, string& b)
{
	//fprintf(stderr, "addSub(%s,%s)\n", a.c_str(), b.c_str());
	string ac = reduce(a);
	string bc = reduce(b);
	if (ac == bc)
		return;
	int ap = ac.length();
	int bp = bc.length();
	if (ap < bp) {
		substitutions.push_back(Substitution(bc, ac));
	} else if (ap > bp) {
		substitutions.push_back(Substitution(ac, bc));
	} else {
    int ax = x_pow(ac);
    int bx = x_pow(bc);
		if (ax < bx) {
		  substitutions.push_back(Substitution(bc, ac));
    } else if (ax > bx) {
		  substitutions.push_back(Substitution(ac, bc));
    } else {
      int ay = y_pow(ac);
      int by = y_pow(bc);
      if (ay < by) {
        substitutions.push_back(Substitution(bc, ac));
      } else if (ay > by) {
        substitutions.push_back(Substitution(ac, bc));
      } 
    }
	}
//	printf("added %s -> %s\n", substitutions.back().s.c_str(),
//		substitutions.back().rep.c_str());
  return;
}

void CanonicalName::add_relator(string relator)
{
	relators.push_back(relator);
}

void CanonicalName::add_relator_internal(string relator)
{
	//fprintf(stderr, "add_relator_internal(%s)\n", relator.c_str());
	string rr(relator + relator);
	string::size_type l = relator.size();
	string::size_type sl;
	string::size_type pos;
	for (pos = 0; pos < l; ++pos) {
		for (sl = 1; sl < l; ++sl) {
			string a(rr.substr(pos, sl));
			string b(inverse(rr.substr(pos+sl, l-sl)));
			add_substitution(a, b);
		}
	}
}

#define MAX_LEN 200

string CanonicalName::reduce(string s)
{
	set<string> visited;
  string ms(s);
	bool done = false;
	while (!done) {
		done = true;
		vector<Substitution>::iterator it;
		for (it = substitutions.begin(); it != substitutions.end(); ++it) {
			string::size_type pos = s.find(it->s);
			if (pos != string::npos) {
				done = false;
				// fprintf(stderr, "replacing in %s: %s -> %s\n", s.c_str(), it->s.c_str(), it->rep.c_str());
				s.replace(pos, it->s.length(), it->rep);
        if (s.length() < ms.length()) {
          ms = s;
        }
				if (visited.find(s) != visited.end() || s.length() > MAX_LEN) {
					fprintf(stderr, "loop detected in canonical name reduce %s to %s\n", ms.c_str(), s.c_str());
					return ms;
				}
				visited.insert(s);
			}
		}
	}
	return s;
}

string CanonicalName::get_canonical_name(string s)
{
	if (!impl) {
		initImpl();
	}
	string result = impl->reduce(inverse(impl->reduce(inverse(s))));
	//fprintf(stderr, "CanonicalName(%s) = %s\n", s.c_str(), result.c_str());
	return result;
}
