#include "PartialTree.hh"
#include "TestCollection.hh"

extern TestCollection g_tests;

// Consume tree from stdin. The tree must be
// provided in pre-order depth-first traversal.
PartialTree read_tree()
{
  PartialTree t;
  char buf[1000];
  if (!fgets(buf, sizeof(buf), stdin)) {
    fprintf(stderr, "unexpected EOF\n");
    exit(1);
  }
  int n = strlen(buf);
  if (buf[n-1] == '\n')
    buf[n-1] = '\0';
  if (buf[0] == 'X') {
    t.result.index = -1;
    t.l_child = new PartialTree(read_tree());
    t.r_child = new PartialTree(read_tree());
  } else if (strstr(buf, "HOLE") != NULL) {
    t.result.index = -2;
  } else {
    if (isdigit(buf[0])) {
      t.result.index = atoi(buf);
    } else {
      if (strchr("xXyY", buf[2]) != NULL) {
        t.result.index = g_tests.add(std::string(buf));
      }
    }
    if (g_debug) {
      fprintf(stderr, "%d\n", t.result.index);
    }
  }
  return t;
}

void truncate_tree(PartialTree& t)
{
  if (t.l_child) {
    truncate_tree(*t.l_child);
    delete t.l_child;
    t.l_child = 0;
  }
  if (t.r_child) {
    truncate_tree(*t.r_child);
    delete t.r_child;
    t.r_child = 0;
  }
}

int tree_size(PartialTree& t) {
  int size = 1;
  if (t.l_child)
    size += tree_size(*t.l_child);
  if (t.r_child)
    size += tree_size(*t.r_child);
  return size;
}
