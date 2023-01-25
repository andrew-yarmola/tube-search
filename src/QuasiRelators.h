#ifndef __QUASI_RELATORS_H
#define __QUASI_RELATORS_H

#include <string>
#include <map>
#include <set>
#include <vector>
#include "types.hh"
#include "SL2.hh"
#include "CanonicalName.hh"

class QuasiRelators {
public:
	std::string get_name(std::string w);               // get the canonical name of a quasi-relator
	void add_quasi_relator(std::string w);        // record that this word is a quasi-relator

	std::vector<std::string> all_words();
	std::vector<std::string> word_classes();
	std::string desc();                               // string describing this set of quasi-relators
	std::string min_pow_desc();                               // string describing minimal power set of quasi-relators
	bool is_quasi_relator(std::string w);         // is this word a quasi-relator?

  template<typename T>
  std::string desc(const Params<T>& p);

private:
  CanonicalName canonical_name;
	typedef std::map< std::string, std::string > NameStore;
	NameStore names;
	std::vector<std::string> name_vector;
	std::string inverse(std::string w);
};

template<typename T>
std::string QuasiRelators::desc(const Params<T>& p)
{
  std::string buf;
  std::string word;
  std::set<std::string> words;
	for (auto qr : name_vector) {
		if (!buf.empty() && buf.back() != ',')
			buf += ",";
    word = proven_identity(qr, p);
    word = canonical_name.get_canonical_name(word);
    if (word.length() > 0) {
      if (words.insert(word).second) {
        buf += word;
      }
    }
	}
	return buf;
}

#endif
