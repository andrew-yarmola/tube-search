#ifndef __partialtree_h
#define __partialtree_h

#include "types.hh"

struct PartialTree {
  PartialTree *l_child = NULL;
  PartialTree *r_child = NULL;
  TestResult result = {-1, open, word_pair()};
};

// Consume tree from stdin. The tree must be
// provided in pre-order depth-first traversal.
PartialTree read_tree();

void truncate_tree(PartialTree& t);

int tree_size(PartialTree& t);

#endif // __partialtree_h
