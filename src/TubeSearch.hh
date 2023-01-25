#include <vector>
#include <set>
#include <utility>
#include <string>
#include "types.hh"

std::vector< word_pair > find_pairs(
	Params<Complex> center,
	std::vector< std::string > seed_words,
	int num_words,
	int max_length,
	std::vector< std::string > quasi_relators
);

std::vector< word_pair > find_words_v2(
          const Params<Complex>& params,
          int num_words, int max_move_len,
          const std::vector<std::string>& relators, const std::map<std::string, int>& seen);
