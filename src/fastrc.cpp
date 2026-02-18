#include <Rcpp.h>
using namespace Rcpp;

static unsigned char dna_comp_table[256];
static bool initialized = false;

void init_table(bool rna_mode) {
  for (int i = 0; i < 256; i++)
    dna_comp_table[i] = (unsigned char)i;

  // Standard mappings
  dna_comp_table['C'] = 'G';
  dna_comp_table['c'] = 'g';
  dna_comp_table['G'] = 'C';
  dna_comp_table['g'] = 'c';
  dna_comp_table['U'] = 'A';
  dna_comp_table['u'] = 'a';
  dna_comp_table['T'] = 'A';
  dna_comp_table['t'] = 'a';

  // Conditional mapping for A
  if (rna_mode) {
    dna_comp_table['A'] = 'U';
    dna_comp_table['a'] = 'u';
  } else {
    dna_comp_table['A'] = 'T';
    dna_comp_table['a'] = 't';
  }

  // IUPAC Ambiguities
  dna_comp_table['M'] = 'K';
  dna_comp_table['m'] = 'k';
  dna_comp_table['K'] = 'M';
  dna_comp_table['k'] = 'm';
  dna_comp_table['R'] = 'Y';
  dna_comp_table['r'] = 'y';
  dna_comp_table['Y'] = 'R';
  dna_comp_table['y'] = 'r';
  dna_comp_table['S'] = 'S';
  dna_comp_table['s'] = 's';
  dna_comp_table['W'] = 'W';
  dna_comp_table['w'] = 'w';
  dna_comp_table['V'] = 'B';
  dna_comp_table['v'] = 'b';
  dna_comp_table['B'] = 'V';
  dna_comp_table['b'] = 'v';
  dna_comp_table['H'] = 'D';
  dna_comp_table['h'] = 'd';
  dna_comp_table['D'] = 'H';
  dna_comp_table['d'] = 'h';
  dna_comp_table['N'] = 'N';
  dna_comp_table['n'] = 'n';

  initialized = true;
}

// [[Rcpp::export]]
CharacterVector fast_rc_impl(CharacterVector seqs, std::string type = "DNA") {
  bool rna_mode = (type == "RNA");

  // Re-initialize only if mode changes or first run
  static std::string current_type = "";
  if (!initialized || current_type != type) {
    init_table(rna_mode);
    current_type = type;
  }

  int n = seqs.size();
  CharacterVector res(n);

  for (int i = 0; i < n; i++) {
    if (CharacterVector::is_na(seqs[i])) {
      res[i] = NA_STRING;
      continue;
    }

    const char *s = seqs[i];
    int len = strlen(s);
    std::string rc_s(len, ' ');

    for (int j = 0; j < len; j++) {
      rc_s[j] = dna_comp_table[(unsigned char)s[len - 1 - j]];
    }
    res[i] = rc_s;
  }
  return res;
}
