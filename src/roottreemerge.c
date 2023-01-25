#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define BUFSZ 1024

using namespace std;

bool g_verbose = false;

// print smallest merged tree
void print_merge(FILE** trees, vector<bool> states, int num_trees, bool do_print)
{
  bool all_done = true;
  bool do_print_after = do_print;
  char buf[BUFSZ];
  char cond[BUFSZ];
  cond[0] = '\0';
  for (int i = 0; i < num_trees; ++i) {
    if (states[i]) { 
      if (!fgets(buf, sizeof(buf), trees[i])) {
        fprintf(stderr, "unexpected EOF\n");
        exit(1);
      }
      if (buf[0] != 'X') {
        states[i] = false;
        if (buf[0] != '6') {
          strncpy(cond, buf, BUFSZ);
          do_print_after = false;
        }
      } else {
        all_done = false;
      }
    }
  }
  if (!all_done) {
    if (do_print_after) {
      printf("X\n");
    } else if (do_print) {
      printf(cond);
    }
    // should deep copy the states vector
    print_merge(trees, states, num_trees, do_print_after); // left 
    print_merge(trees, states, num_trees, do_print_after); // right
  } else if (do_print) {
    if (strlen(cond) > 0) {
      printf(cond);
    } else {
      printf(buf);
    }
  }
}

FILE* open_box_file(const char* boxcode, char* tree_location)
{
  if (g_verbose) fprintf(stderr, "open box %s at %s\n", boxcode, tree_location);
  FILE* fp;
  char fileboxcode[BUFSZ];
  char file_name[BUFSZ];
  strncpy(fileboxcode, boxcode, BUFSZ);
  if (strncmp(fileboxcode, "", BUFSZ) == 0) {
    strncpy(fileboxcode, "root", BUFSZ);
  }
  snprintf(file_name, BUFSZ, "%s/%s.out", tree_location, fileboxcode);
  struct stat sb;
  if (0 == stat(file_name, &sb)) {
    if (g_verbose) fprintf(stderr, "opening %s\n", file_name);
    fp = fopen(file_name, "r");
    return fp;
  }
  snprintf(file_name, BUFSZ, "%s/%s.out.tar.gz", tree_location, fileboxcode);
  if (0 == stat(file_name, &sb)) {
    char command_buf[BUFSZ];
    if (g_verbose) fprintf(stderr, "opening %s\n", file_name);
    snprintf(command_buf, BUFSZ, "gzcat %s", file_name);
    fp = popen(command_buf, "r");
    return fp;
  }
  return NULL;
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    fprintf(stderr, "Usage: roottreemerge tree_locations(s)\n");
    exit(1);
  }

  int num_trees = argc - 1;
  char** locations = argv + 1;
  FILE** trees = (FILE**) calloc(num_trees, sizeof(FILE*));
  vector<bool> states;
  for (int i = 0; i < num_trees; ++i) {
    FILE* tree = open_box_file("", locations[i]);
    if (tree) {
      trees[i] = open_box_file("", locations[i]);
      states.push_back(true);
    } else {
      fprintf(stderr, "No root tree at %s\n", locations[i]);
      exit(1);
    }
  }
  print_merge(trees, states, num_trees, true);
  // no need to free mem at end
}
