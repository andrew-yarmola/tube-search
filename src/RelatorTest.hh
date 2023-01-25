/*
 *  RelatorTest.h
 *  mom
 *
 *  Created by Nathaniel Thurston on 13/10/2007.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include <string>
#include <vector>
#include <set>
#include "CanonicalName.hh"

struct RelatorTest
{
public:
	// if the return value is false, sets mandatory to
  // the list of subwords which should not be identies.
	bool is_impossible(std::string word,
      std::vector<std::string>& required_non_identities);
	bool is_good(std::string word);
	static RelatorTest* create(const char* impos_path,
                             const char* bad_rel_path);
private:
  void load(const char* impos_path,
            const char* bad_rel_path);
  CanonicalName canonical_name;
  std::set<std::string> always_impossible;
  std::set<std::string> bad_relators; 
};
