#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define BUFSZ 1024
#define PARSE_DEPTH_MAX 12
#define PARSE_SIZE_MAX 100000

bool g_verbose = false;

struct partial_tree {
  partial_tree() : left_child(NULL), right_child(NULL), data() {}
  partial_tree *left_child;
  partial_tree *right_child;
  std::string data;
};

void free_tree(partial_tree* t)
{
  if (t) {
    free_tree(t->left_child);
    free_tree(t->right_child);
  delete t;
  }
}

void truncate_tree(partial_tree& t)
{
  t.data.assign("H");
  if (t.left_child) {
    free_tree(t.left_child);
    t.left_child = NULL;
  }
  if (t.right_child) {
    free_tree(t.right_child);
    t.right_child = NULL;
  }
}

FILE* read_to_box(const char* boxcode, char* tree_location);

// Consume tree from a file stream. The tree must be
// provided in pre-order depth-first traversal.
int read_tree(FILE* fp, std::string& box, partial_tree& t,
    int depth, int target_size, int truncate_depth, char* tree_location)
{
  int subtree_size = 0; // counts leaf nodes
  char buf[BUFSZ];
  if (!fgets(buf, sizeof(buf), fp)) {
    fprintf(stderr, "unexpected EOF\n");
    exit(1);
  }
  if (buf[0] == 'X') {
    t.left_child = new partial_tree();
    t.right_child = new partial_tree();
    box.push_back('0');
    subtree_size += read_tree(fp, box, *t.left_child, depth+1,
      target_size, truncate_depth, tree_location);
    box.pop_back();
    box.push_back('1');
    subtree_size += read_tree(fp, box, *t.right_child, depth+1,
      target_size, truncate_depth, tree_location);
    box.pop_back();
  } else if (buf[0] == 'H') {
    FILE* hole_head = read_to_box(box.c_str(), tree_location);
    if (!hole_head) {
      t.data.assign(buf);
      return 1;
    } else {
      subtree_size += read_tree(hole_head, box, t, depth,
          target_size, truncate_depth, tree_location);
      fclose(hole_head);
    } 
  } else {
    t.data.assign(buf);
    return 1;
  }
  if (subtree_size > target_size  && depth > truncate_depth) {
    truncate_tree(t);  
  }
  // reports all seen, not the truncated size
  return subtree_size;
}

void print_merge(partial_tree** trees, int num_trees)
{
  bool all_holes = true;
  bool all_done = true;
  char data[BUFSZ];
  strncpy(data, "X_SPLIT", BUFSZ);
  for (int i = 0; i < num_trees; ++i) {
    partial_tree* t = trees[i];
    if (t) {
      if (t->data.length() > 0) {
        if (t->data[0] != 'H') {
          all_holes = false;
          if (t->data[0] != '6') {
            strncpy(data, t->data.c_str(), BUFSZ);
            all_done = true;
            break;
          }
        }
      } else {
        all_done = false;
      }
    }
  }
  if (all_done) {
    if (all_holes) {
      printf("HOLE\n");
    } else if (data[0] == 'X') {
      printf("6\n");
    } else {
      printf("%s", data);
    }
  } else if (data[0] == 'X') {
    printf("X\n");
    partial_tree** all_left = (partial_tree**) calloc(num_trees, sizeof(partial_tree*));
    for (int i = 0; i < num_trees; ++i) {
      if (trees[i]) {
        all_left[i] = trees[i]->left_child;
      } else {
        all_left[i] = NULL;
      }
    }
    print_merge(all_left, num_trees);
    free(all_left); 
    partial_tree** all_right = (partial_tree**) calloc(num_trees, sizeof(partial_tree*));
    for (int i = 0; i < num_trees; ++i) {
      if (trees[i]) {
        all_right[i] = trees[i]->right_child;
      } else {
        all_right[i] = NULL;
      }
    }
    print_merge(all_right, num_trees);
    free(all_right); 
  }
}

FILE* open_box_file(char* boxcode, char* tree_location)
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

FILE* read_to_box(const char* boxcode, char* tree_location)
{
  char fileboxcode[BUFSZ];
  strncpy(fileboxcode, boxcode, BUFSZ);

  FILE* fp = 0;
  size_t code_len = strlen(fileboxcode);
  while (code_len >= 0) {
    fileboxcode[code_len] = '\0';
    fp = open_box_file(fileboxcode, tree_location);
    if (fp) { break; }
    --code_len;
  }
  if (!fp) { 
    fprintf(stderr, "Failed to find box %s at %s\n", boxcode, tree_location);
    exit(2);
  }
  char buf[BUFSZ];
  char code_tail[BUFSZ];
  std::string box(fileboxcode); 
  strncpy(code_tail, boxcode + code_len, BUFSZ);
  char* code = code_tail;
  while (*code && fgets(buf, sizeof(buf), fp)) {
    if (buf[0] != 'X') { // If not a split, return a bad FILE
      fclose(fp);
      return NULL;
    }
    box.push_back(*code);
    // Actually have to process the tree if we go right at any point in the boxcode
    if (*code == '1') {
      partial_tree* right = new partial_tree();
      read_tree(fp, box, *right, 0, PARSE_SIZE_MAX, PARSE_DEPTH_MAX, tree_location);
      free_tree(right); // just need the traversal 
    }
    ++code;  // Keeps going left in the tree as *boxcode == 0
  }
  return fp;
}

int main(int argc, char** argv)
{
  if (argc < 5) {
    fprintf(stderr, "Usage: treemerge size truncate_depth boxcode tree_locations(s)\n");
    exit(1);
  }

  int target_size = atoi(argv[1]);
  int truncate_depth = atoi(argv[2]);
  char boxcode[BUFSZ];
  strncpy(boxcode, argv[3], BUFSZ);

  int num_trees = argc - 4;
  char** locations = argv + 4;
  partial_tree** trees = (partial_tree**) calloc(num_trees, sizeof(partial_tree*));    
  for (int i = 0; i < num_trees; ++i) {
    FILE* tree_head = read_to_box(boxcode, locations[i]);
    if (!tree_head) {
      trees[i] = NULL;
    } else {
      partial_tree* t = new partial_tree();
      std::string box(boxcode);
      read_tree(tree_head, box, *t, 0, target_size, truncate_depth, locations[i]);
      trees[i] = t;
    }
  }
  print_merge(trees, num_trees);
  // no need to free mem at end
}
