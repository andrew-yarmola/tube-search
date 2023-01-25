#ifndef __refine_h
#define __refine_h

#include "Box.h"
#include "PartialTree.hh"

struct Options {
  Options() :
  box_name(""), // Binary representation of box
  words_file("words"), // Previously generated words
  impossible_file("null"), // impossible relators 
  bad_relator_file(
      "/u/yarmola/margsym/margulis-search/bad_relators"), // relators that can't be used 
  max_depth(24), // Maximum depth for a file
  truncate_depth(6), 
  invent_depth(12),
  max_size(1000000),
  improve_tree(false),
  word_search_depth(-1),
  fill_holes(false),
  max_word_length(40) {}
  const char* box_name;
  const char* words_file;
  const char* impossible_file;
  const char* bad_relator_file;
  int max_depth;
  int truncate_depth;
  int invent_depth;
  int max_size;
  bool improve_tree;
  int word_search_depth;
  bool fill_holes;
  int max_word_length;
};

void refine_tree(Box box, PartialTree& t);

void print_tree(PartialTree& t);

#endif // __refine_h
