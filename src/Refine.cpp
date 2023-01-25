#include "Refine.hh"
#include "TestCollection.hh"
#include "TubeSearch.hh"
#include "QuasiRelators.h"

using namespace std;

typedef vector< vector< box_state > > TestHistory;

Options g_options;
TestCollection g_tests;
int g_boxes_visited = 0;

#define IMPROVE_MOD 4 
#define IMPROVE_HIST 7 
#define QR_MOD 4
#define WORD_SEARCH_MOD 6
#define WORD_SEARCH_DEPTH 50

extern double g_cosh_marg_upper_bound;
extern double g_cosh_marg_lower_bound;
extern double g_sinh_d_bound; 

extern int num_bound_tests;

extern bool g_debug;

unordered_map<string, SL2<AJ> > short_words_cache;

bool refine_recursive(Box box, PartialTree& t, int depth,
    TestHistory& history, vector< Box >& place,
    int new_depth, int& searched_depth)
{
  place.push_back(box);

  if (t.result.index == -2 && !g_options.fill_holes) {
    return true;
  }

  int old_result_index = t.result.index;
  if (t.result.index >= 0) {
    t.result = g_tests.evaluate_box(t.result.index, box);
    if (t.result.state != open && t.result.state != open_with_qr) {
      return true;
    } else { 
      fprintf(stderr,
          "FAILED to eliminate %s with test %s with result %d\n",
          box.name.c_str(), g_tests.get_name(old_result_index).c_str(),
          t.result.state);
    }
  }
  
  // Check if the box is now small enough that some former qrs actually kill it
  if (depth % QR_MOD == 1) {
    for (auto qr : box.qr.word_classes()) {
      word_pair qr_pair(qr, ""); 
      t.result = g_tests.evaluate_AJ(qr_pair, box);
      if (t.result.state != open && t.result.state != open_with_qr) {
        return true;
      }
    }
  }

  if (g_options.improve_tree || !t.l_child) {
    for (int i = 0; i < g_tests.size(); ++i) {
      // only do boundary tests on a regular basis
      if (i >= num_bound_tests && depth % IMPROVE_MOD == 0) {
        break;
      } 
      vector<box_state>& th = history[i];
      while (th.size() <= depth) {
        box_state result = g_tests.evaluate_center(i, place[th.size()]);
        th.push_back(result);
      }
      bool do_eval = true;
      int s = th.size();
      for (int j = 1; j <= min(s, IMPROVE_HIST); j++) {
        do_eval = do_eval && th[th.size() - j] != open;
      }
      if (do_eval) {
        t.result = g_tests.evaluate_box(i, box);
        if (t.result.state != open && t.result.state != open_with_qr) {
          return true;
        }
      }
    }
  }

  // WORD SEARCH IS ON
  if (g_options.word_search_depth > 0 && depth > 0
      && (g_options.improve_tree || !t.l_child)
      && box.name.length() > g_options.word_search_depth
      && depth % WORD_SEARCH_MOD == 0) {
    Box& search_place = box;
    vector<word_pair> search_pairs_v1 = find_pairs(search_place.center(),
        vector<string>(), 1, g_options.max_word_length, box.qr.word_classes());
    vector<word_pair> search_pairs_v2 = find_words_v2(search_place.center(),
        1, 14, box.qr.word_classes(), g_tests.seen_words);
    vector<word_pair> search_pairs;
    search_pairs.insert(search_pairs.end(),
        search_pairs_v1.begin(), search_pairs_v1.end());
    search_pairs.insert(search_pairs.end(),
        search_pairs_v2.begin(), search_pairs_v2.end());
    while (search_pairs.size() > 0) {
      word_pair new_pair = search_pairs.back();
      int old_size = g_tests.size();
      int new_index = g_tests.add(new_pair);
      history.resize(g_tests.size());
      search_pairs.pop_back();
      if (old_size < g_tests.size()) {
        fprintf(stderr, "search (%s) found (%s,%s) at (%s)\n",
                search_place.qr.desc(box.cover()).c_str(), new_pair.first.c_str(),
                new_pair.second.c_str(), search_place.name.c_str());
        t.result = g_tests.evaluate_box(new_index, box);
        if (t.result.state != open && t.result.state != open_with_qr) {
          return true;
        }
      }
    }
  }

  if (box.qr.word_classes().size() > 0) {
    t.result.state = open_with_qr;
  } else {
    t.result.state = open;
  }

  if (!t.l_child) {
    if (depth >= g_options.max_depth || ++g_boxes_visited >= g_options.max_size
        || ++new_depth > g_options.invent_depth) {
      fprintf(stderr, "HOLE %s (%s)\n", box.name.c_str(),
          box.qr.desc(box.cover()).c_str());
      return false;
    }
    t.l_child = new PartialTree();
    t.r_child = new PartialTree();
  }

  bool is_complete = true;

  is_complete = refine_recursive(box.child(0), *t.l_child, depth + 1, history,
      place, new_depth, searched_depth) && is_complete;
  if (place.size() > depth + 1)
    place.resize(depth + 1);
  for (int i = 0; i < g_tests.size(); ++i) {
    if (history[i].size() > depth)
      history[i].resize(depth);
  }
  if (searched_depth > depth)
    searched_depth = depth;
  if (is_complete || depth < g_options.truncate_depth)
    is_complete = refine_recursive(box.child(1), *t.r_child, depth + 1, history,
        place, new_depth, searched_depth) && is_complete;
  if (old_result_index >= 0 && t.result.index != old_result_index) {
    fprintf(stderr, "invalid box %s(%s) %d %s\n",
        g_tests.get_name(old_result_index).c_str(), box.name.c_str(),
      tree_size(t), is_complete ? "Patched" : "Unpatched");
  }
  if (!is_complete && depth >= g_options.truncate_depth) {
    truncate_tree(t);
  }
  return is_complete;
}

void refine_tree(Box box, PartialTree& t)
{
  TestHistory history(g_tests.size());
  vector<Box> place;
  int searched_depth = 0;
  refine_recursive(box, t, 0, history, place, 0, searched_depth);
}

void print_tree(PartialTree& t)
{
    char type = 'F';
    word_pair p = t.result.words;
    switch (t.result.state) {
      case open :
      case open_with_qr : {
        if (t.l_child && t.r_child) {
          printf("X\n");
          print_tree(*t.l_child);
          print_tree(*t.r_child);
        } else {
          printf("HOLE\n");
        }
        return;
      }
      case killed_bounds : {
        printf("%s\n", g_tests.get_name(t.result.index).c_str());
        return;
      }
      case killed_impossible_relator : type = 'E'; break; 
      case killed_x_hits_y : type = 'a'; break;
      case killed_y_hits_x : type = 'A'; break;
      case killed_x_hits_x : type = 'x'; break;
      case killed_y_hits_y : type = 'y'; break;
      case killed_x_not_cyclic : type = 'p'; break;
      case killed_y_not_cyclic : type = 'P'; break;
      case killed_move : type = 'm'; break;
      case killed_marg : type = 'M'; break;
      case killed_nbd_x : type = 'c'; break;
      case killed_nbd_y : type = 'C'; break;
      case killed_nbd : type = 'D'; break;
      case killed_w_ax_hits_sym_axis: type = 'n'; break;
      case killed_w_ay_hits_sym_axis: type = 'N'; break;
      case killed_sym: type = 'S'; break; 
      case proven_relator: type = 'R'; break; 
      case var_nbd_x : type = 'v'; break;
      case var_nbd_y : type = 'V'; break;
      case var_nbd : type = 'W'; break;
      default : return;
    }
    
    printf("%c(%s,%s)\n", type, p.first.c_str(), p.second.c_str()); 
}

