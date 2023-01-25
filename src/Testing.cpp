#include "Refine.hh"
#include "TestCollection.hh"
#include "TubeSearch.hh"
#include "QuasiRelators.h"

using namespace std;

typedef vector< vector< box_state > > TestHistory;

// for TestCollection extern compatibility
double g_cosh_marg_upper_bound = 1.2947;
double g_cosh_marg_lower_bound = 1.0054;
double g_sinh_d_bound = 1.3426; 

Options g_options;
int g_boxesVisited = 0;

TestCollection g_tests;

box_state evaluate_center(TestCollection& tests, int index, Box& box)
{
  //  fprintf(stderr, "Evaluating box test index %d\n", index);
  Params<Complex> center = box.center();
  switch(index) {
    case 0:	{
              if (absLB(abs(center.coshmu - 1) + 
                  abs(sqrt(center.sintx2 + 1.1)) +
                  abs(sqrt(center.sinty2) - 0.5) +
                  abs(center.sinhdy + 0.7)) > 0 || 
                  absLB(abs((center.cosf - 1.2) * (center.sinhdx - 0.5)) +
                  abs((center.sinhdx - 0.3) / (center.sinhdy + 0.6))) > 0) {
                return killed_bounds;
              }
                return open;
            } 
    default: {
               return  open;
             }
  }
}

box_state evaluate_box(TestCollection& tests, int index, Box& box)
{
  //  fprintf(stderr, "Evaluating box test index %d\n", index);
  Params<AJ> cover = box.cover();
  switch(index) {
    case 0:	{
              if (absLB(abs(cover.coshmu - 1) + 
                  abs(sqrt(cover.sintx2 + 1.1)) +
                  abs(sqrt(cover.sinty2) - 0.5) +
                  abs(cover.sinhdy + 0.7)) > 0 || 
                  absLB(abs((cover.cosf - 1.2) * (cover.sinhdx - 0.5)) +
                  abs((cover.sinhdx - 0.3) / (cover.sinhdy + 0.6))) > 0) {
                return killed_bounds;
              }
                return open;
            } 
    default: {
               return  open;
             }
  }
}


bool refine_recursive(Box box, PartialTree& t, int depth, TestHistory& history, vector< Box >& place, int newDepth, int& searched_depth)
{
  //fprintf(stderr, "rr: %s depth %d placeSize %lu\n", box.name.c_str(), depth, place.size());
  place.push_back(box);
  int old_test_index = t.test_index;
  vector<string> new_qrs;

  if (t.test_index >= 0) {
    //    fprintf(stderr, "********************************* Validation *********************************\n");
    //    fprintf(stderr, "%s", box.desc().c_str());
    box_state result = evaluate_box(g_tests, t.test_index, box);
    if (result != open)  {
      t.test_result = result;
      //      fprintf(stderr, "Eliminated %s with test %s with result %d\n", box.name.c_str(), g_tests.get_name(t.test_index), result);
      //      fprintf(stderr, "********************************* End Validation *********************************\n");
      //      fprintf(stderr, "Test %s kills\n %s", g_tests.get_name(t.test_index).c_str(), box.desc().c_str());
      return true;
    } else { 
      fprintf(stderr, "FAILED to eliminate %s with test %s with result %d\n", box.name.c_str(), g_tests.get_name(t.test_index).c_str(), result);
    }
    //    fprintf(stderr, "********************************* End Validation *********************************\n");
  }

  if (t.test_index == -2 && !g_options.fill_holes) {
    return true;
  }

  if (g_options.improve_tree || !t.l_child) {
    for (int i = 0; i < g_tests.size(); ++i) {
      vector<box_state>& th = history[i];
      while (th.size() <= depth && (th.size() < depth-6 || th.empty() || th.back() == open)) {
        //        fprintf(stderr, "********************************* Center Test *********************************\n");
        //        fprintf(stderr, "%s", box.desc().c_str());
        box_state result = evaluate_center(g_tests, i, place[th.size()]);
        //        fprintf(stderr, "********************************* End Center Test *********************************\n");
        th.push_back(result);
      }
      if (th.back() != open) {
        //        fprintf(stderr, "********************************* Evaluate *********************************\n");
        //        fprintf(stderr, "%s", box.desc().c_str());
        box_state result = evaluate_box(g_tests, i, box);
        //        fprintf(stderr, "********************************* End Evaluate *********************************\n");

        switch (result) {
          case killed_bounds: {
                                t.test_index = i;
                                t.test_result = result;
                                return true;
                              }
          default : {
                      continue;
                    }
        }
      }
    }
  }

  t.test_index = -1;

  if (!t.l_child) {
    if (depth >= g_options.max_depth || ++g_boxesVisited >= g_options.max_size || ++newDepth > g_options.invent_depth) {
      //    fprintf(stderr,"Deph %d, max depth %d, boxes_visited %d, max size %d, newDepth %d, invent depth %d\n", depth, g_options.max_depth, g_boxesVisited, g_options.max_size, newDepth, g_options.invent_depth);
      fprintf(stderr, "HOLE %s (%s)\n", box.name.c_str(), box.qr.desc().c_str());
      return false;
    }
    t.l_child = new PartialTree();
    t.r_child = new PartialTree();
  }

  bool is_complete = true;

  is_complete = refine_recursive(box.child(0), *t.l_child, depth+1, history, place, newDepth, searched_depth) && is_complete;
  if (place.size() > depth+1)
    place.resize(depth+1);
  for (int i = 0; i < g_tests.size(); ++i) {
    if (history[i].size() > depth)
      history[i].resize(depth);
  }
  if (searched_depth > depth)
    searched_depth = depth;
  if (is_complete || depth < g_options.truncate_depth)
    is_complete = refine_recursive(box.child(1), *t.r_child, depth+1, history, place, newDepth, searched_depth) && is_complete;
  if (old_test_index >= 0 && t.test_index != old_test_index) {
    fprintf(stderr, "invalid box %s(%s) %d %s\n", g_tests.get_name(old_test_index).c_str(), box.name.c_str(),
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
  //  printf("%d\n", t.test_result);
  switch (t.test_result) {
    case open :
    case open_with_qr : {
                          if (t.l_child && t.r_child) {
                            printf("X\n");
                            print_tree(*t.l_child);
                            print_tree(*t.r_child);
                          } else {
                            printf("HOLE (%s)\n", t.qr_desc.c_str());
                          }
                          return;
                        }
    case killed_bounds : {
                           printf("%s\n", g_tests.get_name(t.test_index).c_str());
                           return;
                         }
    case killed_failed_qr : {
                              printf("%c(%s,)\n", 'Q', t.aux_word.c_str());
                              return;
                            }
    case killed_only_elliptic : type = 'E'; break; 
    case killed_x_hits_y : type = 'a'; break;
    case killed_y_hits_x : type = 'A'; break;
    case killed_x_tube : type = 'x'; break;
    case killed_y_tube : type = 'y'; break;
    case killed_lox_not_x_power : type = 'p'; break; 
    case killed_lox_not_y_power : type = 'P'; break;
    case killed_move : type = 'm'; break;
    case killed_marg : type = 'M'; break;
    case variety_nbd_x : type = 'v'; break;
    case variety_nbd_y : type = 'V'; break;
    case variety_nbd : type = 'W'; break;
    case var_x_hits_y : type = 'c'; break;
    case var_y_hits_x : type = 'C'; break;
    default : return;
  }
  printf("%c%s\n", type, g_tests.get_name(t.test_index).c_str());
}

