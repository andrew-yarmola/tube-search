#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <dirent.h>
#include <map>
#include <vector>
#include <string>
#include <getopt.h>
#include <inttypes.h>

using namespace std;

struct Config {
  bool print_killed = false;
  bool recursive = false;
  bool verbose = false;
  bool print_stats = true;
  bool start_is_root = false;
  char* tree_location;
  char kill_test[1000];
};

#define DIM 6
#define SCL 8   
#define BAL 16 
double bounding_volume;

struct Config g_config;

vector<string> unopened_out_files;
map<string, uint64_t> elimination_counts;
map<string, double> elimination_volumes;
map<string, uint64_t> type_counts;
map<string, double> type_volumes;
uint64_t terminal_node_count = 0;
uint64_t internal_nodes = 0;
double total_volume = 0;

FILE* open_box(char* boxcode, char* file_name)
{
	FILE* fp;
	char boxcode_file[10000];
	strcpy(boxcode_file, boxcode);
  // Open the root file if empty
	if (strcmp(boxcode_file, "") == 0) {
		strcpy(boxcode_file, "root");
  }
	sprintf(file_name, "%s/%s.out", g_config.tree_location, boxcode_file);
	struct stat sb;
	if (0 == stat(file_name, &sb)) {
		if (g_config.verbose) {
      fprintf(stderr, "opening %s\n", file_name);
    }
    fp = fopen(file_name, "r");
    return fp;
	}
  // Look for a gzipped file
	sprintf(file_name, "%s/%s.out.tar.gz", g_config.tree_location, boxcode_file);
	if (0 == stat(file_name, &sb)) {
		char command_buf[10000];
		if (g_config.verbose) {
      fprintf(stderr, "opening %s\n", file_name);
    }
		sprintf(command_buf, "tar -xOzf %s", file_name); // cool trick
		fp = popen(command_buf, "r");
		return fp;
	}
	return 0;
}

bool ends_with(const char *str, const char *suffix)
{
  if (!str || !suffix) {
    return false;
  }
  size_t lenstr = strlen(str);
  size_t lensuffix = strlen(suffix);
  if (lensuffix > lenstr) {
    return false;
  }
  return strncmp(str + lenstr - lensuffix, suffix, lensuffix) == 0;
}

#define BUF_SIZE 1048576

bool  put_stream(FILE* dest, FILE* source) {
  size_t size;
  char buf[BUFSIZ];
  size_t total = 0;
  while ((size = fread(buf, 1, BUFSIZ, source))) {
    fwrite(buf, 1, size, dest);
    total += size;
  }
  if (ferror(dest) != 0) {
    fprintf(stderr, "failed destination: bytes %d/%d\n", total, BUFSIZ);
    return false;
  } else if (ferror(source) != 0) {
    fprintf(stderr, "failed source: bytes %d/%d\n", total, BUFSIZ);
    return false;
  } else {
    return true;
  }
}

bool process_tree(FILE* fp, FILE* out, char* boxcode, double vol) {
	int boxdepth = strlen(boxcode);
	char buf[10000];
	char tmp[10000];
  char file_name[10000];
	int depth = 0;
	while (fgets(buf, sizeof(buf), fp)) {
		bool hole_filled = false;
		if (buf[0] == 'H') {
		  if (g_config.recursive && depth > 0) {
        FILE* fp_hole = open_box(boxcode, file_name);
        if (fp_hole) {
          FILE* hole_out = tmpfile();
          if (hole_out) {
            bool success = process_tree(fp_hole, hole_out, boxcode, vol);
            if (!success) {
              fprintf(stderr, "Error with hole at %s\n", boxcode);
            } else {
              hole_filled = true;
            }
            fclose(fp_hole);
            fclose(hole_out);
          }
        }
      }
		}
    size_t end = strlen(buf);
    strncpy(tmp, buf, end);
    tmp[end-1] = '\0';
    if (g_config.print_stats && !hole_filled) {
      if (buf[0] == 'X') {
        internal_nodes += 1;
      } else {
        string key(tmp);  
        if (elimination_counts.find(key) == elimination_counts.end()) {
          elimination_counts[key] = 0;
          elimination_volumes[key] = 0;
        }
        elimination_counts[key] += 1;
        elimination_volumes[key] += vol; 
        terminal_node_count += 1;
        total_volume += vol; 
      }
    }
    if (g_config.print_killed && 
        strncmp(g_config.kill_test, buf, strlen(g_config.kill_test)) == 0) {
      fprintf(out, "%s: %s\n", tmp, boxcode); // Print the killed boxcode to stdout
    }
		if (buf[0] == 'X') {
      // Descend via left branch
			boxcode[boxdepth + depth] = '0';
			depth += 1;
      vol *= 0.5;
			boxcode[boxdepth + depth] = '\0';
		} else {
      // Go up as many nodes as necessary
			for (; depth > 0 && boxcode[boxdepth + depth-1] == '1'; --depth) {
        vol *= 2;
      }
			if (depth > 0) {
				boxcode[boxdepth + depth-1] = '1'; // Jump from left to right node
				boxcode[boxdepth + depth] = '\0'; // Truncate to keep box current
			} else {
				boxcode[boxdepth] = '\0'; // Truncate to keep box current
				return true;
			}
		}
	}
    
  // If we get to this point, the tree is incomplete
  fprintf(stderr, "The tree is incomplete, run treecat\n");
  return false; 
}

void usage() {
    fprintf(stderr,
        "Usage: treestat [--killed_by] <test> [-v] [-r] tree_location boxcode\n");
    fprintf(stderr, 
        "killed_by: give all terminal nodes from boxcode killed by the given test\n");
    fprintf(stderr, 
        "v: verbose\n");
    fprintf(stderr, 
        "r: recur over all subtree files to prince full subtree of boxcode\n");
    exit(1);
}
		
static struct option long_options[] = {
  {"recursive", no_argument, NULL, 'r' },
  {"verbose", no_argument, NULL, 'v' },
  {"killed_by", required_argument, NULL, 'k' },
  {NULL, 0, NULL, 0}
};

static char opt_str[1000] = "";
void set_opt_str() {
  char* osp = opt_str;
  for (int i = 0; long_options[i].name; ++i) {
    *osp++ = long_options[i].val;
    if (long_options[i].has_arg != no_argument)
      *osp++ = ':';
  }
  *osp = '\0';
}

int main(int argc, char** argv)
{
  set_opt_str();
  if (argc < 2) {
    usage();
    exit(1);
  }
  int ch;
  while ((ch = getopt_long(argc, argv, opt_str, long_options, NULL)) != -1) {
    switch(ch) {
      case 'r': g_config.recursive = true; break;
      case 'v': g_config.verbose = true; break;
      case 'k': {
        strncpy(g_config.kill_test, optarg, sizeof(g_config.kill_test));
        g_config.print_killed = true;
        g_config.print_stats = false;
        break;
      }
      default: usage(); exit(1);
    }
  }
  double corr = pow(BAL, DIM);
  bounding_volume = pow(2 * SCL * BAL, DIM);
  for (int i = 0; i < DIM; ++i) {
    bounding_volume *= pow(2, -i / float(DIM));
  }
  fprintf(stderr, "Bounding volume for dim %d and scale %d is %f\n",
      DIM, SCL, bounding_volume / corr);

  // The fullboxcode parameter can specify the file_name and sequetial boxcode
  // A boxcode is just a sequence of zeros and ones giving
  // a posiiton in a binary tree depth-first traversal
  // The treeFile will also be in pre-order depth-first
  char fullboxcode[10000];
	char boxcode_file[10000];

  if (optind + 1 < argc) {
	  g_config.tree_location = argv[optind];
    strncpy(fullboxcode,  argv[optind+1], 10000);
    strncpy(boxcode_file, argv[optind+1], 10000);
  } else {
    usage();
    exit(1);
  }
  
	int code_length = strlen(fullboxcode);
  if (code_length == 0 || strncmp(fullboxcode, "root", code_length) == 0) {
    g_config.start_is_root = true;
  }

  // See if a file with the tree for a prefix of the box exists
	FILE* fp = 0;
  char file_name[10000];
	while (code_length >= 0) {
		boxcode_file[code_length] = '\0';
		fp = open_box(boxcode_file, file_name);
		if (fp) {
      break;
    }
    --code_length;
	}

  if (!fp) {
    fprintf(stderr, "no box files at tree location\n");
    exit(1);
  }
    
  char * boxcode_const = (char *)calloc(10000, sizeof(char));
	strncpy(boxcode_const, fullboxcode+code_length, 10000);
  char * boxcode = boxcode_const;

	char buf[10000];
  // If the boxcode is still not empty, we traverse down the tree and print only
  // once we get to the proper node. We terminate early if the node does not exist
	while (*boxcode && fgets(buf, sizeof(buf), fp)) {
		if (buf[0] != 'X') { // If not a splitting, print the test failed by the truncated box
			*boxcode = '\0';
			fprintf(stderr, "terminal box = %s%s\n", boxcode_file, boxcode_const);
      free(boxcode_const);
      fclose(fp);
			exit(0);
		}
		if (*boxcode == '1') { // Actually have to process the tree
      // if we go right at any point in the boxcode
      FILE* dev_null = fopen("/dev/null","w");
      bool true_print_stats = g_config.print_stats;
      bool true_print_killed = g_config.print_killed;
      // turn off for now
      g_config.print_stats = false;
      g_config.print_killed = false;
      int success = process_tree(fp, dev_null, boxcode, 0);
      g_config.print_stats = true_print_stats;
      g_config.print_killed = true_print_killed;
      fclose(dev_null); 
      if (!success) exit(1); // Incomplete tree or boxcode not found
    }
		++boxcode; // Keeps going left in the tree as *boxcode == 0
	}

  free(boxcode_const);
  FILE* out = tmpfile();
  if (!out) {
    exit(1);
  }
  double starting_volume = bounding_volume * pow(0.5, strlen(fullboxcode));
  fprintf(stderr, "Starting volume is %f for box %s\n",
      starting_volume / corr, fullboxcode);
	bool success = process_tree(fp, out, fullboxcode, starting_volume);
  fclose(fp);

  if (!success) {
    fclose(out);
    exit(1);
  } else {
    rewind(out);
    bool put_success = put_stream(stdout, out);
    fclose(out);
    if (!put_success) {
      fprintf(stderr, "failed to write final output\n");
      exit(1);
    }
    if (g_config.print_stats) {
      for (auto& it: elimination_counts) {
        string type = it.first.substr(0,1);
        if (type_counts.find(type) == type_counts.end()) {
          type_counts[type] = 0;
          type_volumes[type] = 0;
        }
        double vol = elimination_volumes[it.first];
        type_counts[type] += it.second;
        type_volumes[type] += vol;
        printf("%s: %" PRIu64 ": %f\n", it.first.c_str(), it.second, vol/corr);
      }
      for (auto& it: type_counts) {
        double vol = type_volumes[it.first];
        printf("%s: %" PRIu64 ": %f\n", it.first.c_str(), it.second, vol/corr);
      }
      printf("internal nodes: %" PRIu64 "\n", internal_nodes);
      printf("total: %" PRIu64 ": %f\n", terminal_node_count, total_volume/corr);
    }
    else exit(0); 
  }
}
