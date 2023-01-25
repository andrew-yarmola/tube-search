/*
 *  RelatorTest.cpp
 *
 *  Created by Nathaniel Thurston on 13/10/2007.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "RelatorTest.hh"
#include "types.hh"
#include <cstdio>

using namespace std;


bool RelatorTest::is_good(string word)
{
  if (word.length() == 0 ||
      bad_relators.find(word) != bad_relators.end()) {
    return false;
  }
  return true;
}

bool RelatorTest::is_impossible(string word,
    vector<string>& required_non_identities)
	{
    required_non_identities.clear();
    string canon = canonical_name.get_canonical_name(word);
    if (syllables(canon) < 5 ||
        always_impossible.find(canon) != always_impossible.end()) {
      return true;
    }
    string cycle(canon);
    int rot = 1;
    // could be optimized
    while (rot != canon.length()) {
      rotate(cycle.begin(), cycle.begin() + 1, cycle.end());
      if (cycle == canon) {
        string sub = canon.substr(0, rot);
        if (syllables(canonical_name.get_canonical_name(sub)) < 5) {
          return true;
        } else {
          required_non_identities.push_back(sub);
		      return false;
        }
      }
      ++rot;
    }
		return false;
	}

void RelatorTest::load(const char* impos_path,
                       const char* bad_rel_path)
{
  char buf[1000];
  char word_buf[1000];
  FILE* fp = fopen(impos_path, "r");
  while (fp && fgets(buf, sizeof(buf), fp)) {
    if (buf[0] != '/') {
      buf[strcspn(buf, "\r\n")] = 0;
      string canon = canonical_name.get_canonical_name(string(buf));
      always_impossible.insert(canon);
    }
  }
  fclose(fp);
  fp = fopen(bad_rel_path, "r");
  while (fp && fgets(buf, sizeof(buf), fp)) {
    if (buf[0] != '/') {
      buf[strcspn(buf, "\r\n")] = 0;
      string canon = canonical_name.get_canonical_name(string(buf));
      bad_relators.insert(canon);
    }
  }
}

RelatorTest* RelatorTest::create(const char* impos_path,
                                 const char* bad_rel_path)
{
	RelatorTest* impossible = new RelatorTest();
	impossible->load(impos_path, bad_rel_path);
	return impossible;
}
