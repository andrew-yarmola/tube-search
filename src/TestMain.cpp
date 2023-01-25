#include <set>
#include <getopt.h>

#include "Testing.hh"
#include "TestCollection.hh"

using namespace std;

extern Options g_options;
extern TestCollection g_tests;
extern int g_boxesVisited;

const char* g_programName;

static struct option long_options[] = {
  {"box",  required_argument, NULL, 'b' },
  {"words", required_argument, NULL, 'w' },
  {"powers", required_argument, NULL, 'p'},
  {"max_depth", required_argument, NULL, 'd' },
  {"invent_depth", required_argument, NULL, 'i' },
  {"improve_tree", no_argument, NULL, 'I'},
  {"truncate_depth", required_argument, NULL, 't' },
  {"max_size", required_argument, NULL, 's' },
  {"word_search_depth", required_argument, NULL, 'B'},
  {"fill_holes", no_argument, NULL, 'f'},
  {"cosh_margulis_upper", required_argument, NULL, 'm'},
  {"sinh_tube_upper", required_argument, NULL, 'r'},
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

void usage()
{
  const char* longUsage = "\
       --box <box_id>\n\
           Box ID for the root of the input and output trees.\n\
    May include characters not in [01], which are ignored.\n\
\n\
Options controlling which relators to use:\n\
  [ --words  <words_file> ]\n\
    File containing starting list of words to try.\n\
  [ --word_search_depth <n> ]\n\
    Perform a search for relators when visiting a node at least n levels deep.\n\
\n\
Options controlling which relators eliminate boxes:\n\
  [ --powers <powers_file> ]\n\
    File containing impossible relator definitions.\n\
    See ImpossibleRelators::load(...)\n\
\n\
Options controlling tree manipulation:\n\
  [ --max_depth <n> ]\n\
    Don't descend more than n levels deeper than the root box.\n\
  [ --invent_depth <n> ]\n\
    Don't descend more than n levels deeper than the terminal node of the input tree.\n\
  [ --max_size <n> ]\n\
    Don't allow the output tree to have more than n nodes.\n\
  [ --truncate_depth <n> ]\n\
    Don't emit holes more than n levels deeper than the root node.\n\
    Instead, replace the subtree-with-holes with a single hole.\n\
  [ --improve_tree ]\n\
    If set, attempt to directly eliminate internal nodes of the input tree.\n\
  [ --fill_holes ]\n\
    If set, attempt to patch holes in the input tree.\n\
";
  fprintf(stderr, "Usage: %s %s\n\n%s", g_programName, opt_str, longUsage);
}

void load_words(set<string>& s, const char* fileName)
{
  FILE* fp = fopen(fileName, "r");
  char buf[10000];
  while (fp && fgets(buf, sizeof(buf), fp)) {
    int n = -1 + strlen(buf);
    if (buf[n] == '\n')
      buf[n] = '\0';
    s.insert(buf);
  }
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
    fprintf(stderr,"Arg %c, %s\n", ch, optarg);
    switch(ch) {
    case 'b': g_options.box_name = optarg; break;
    case 'w': g_options.words_file = optarg; break;
    case 'p': g_options.powers_file = optarg; break;
    case 'd': g_options.max_depth = atoi(optarg); break;
    case 'i': g_options.invent_depth = atoi(optarg); break;
    case 'I': g_options.improve_tree = true; break;
    case 't': g_options.truncate_depth = atoi(optarg); break;
    case 's': g_options.max_size = atoi(optarg); break;
    case 'B': g_options.word_search_depth = atoi(optarg); break;
    case 'f': g_options.fill_holes = true; break;
    }
  }

  Box box;
  for (const char* boxP = g_options.box_name; *boxP; ++boxP) {
    if (*boxP == '0') {
      box = box.child(0);
    } else if (*boxP == '1') {
      box = box.child(1);
    }
  }
 
  g_tests.load(g_options.words_file);
  g_tests.load_impossible_relations(g_options.powers_file);

  fprintf(stderr, "%s", box.desc().c_str());
  PartialTree t = read_tree();
  refine_tree(box, t);
  print_tree(t);
  fprintf(stderr, "%d nodes added\n", g_boxesVisited);
}
