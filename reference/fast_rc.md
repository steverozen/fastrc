# Fast Reverse Complement of DNA/RNA Sequences

Computes the reverse complement of one or more DNA or RNA sequences
using a C++ lookup table for O(1) complement mapping. Supports full
IUPAC ambiguity codes and preserves case.

## Usage

``` r
fast_rc(seqs, type = "DNA")
```

## Arguments

- seqs:

  A character vector of DNA or RNA sequences.

- type:

  Either `"DNA"` (A\\\leftrightarrow\\T) or `"RNA"`
  (A\\\leftrightarrow\\U).

## Value

A character vector of reverse-complemented sequences, with `NA` values
preserved.

## Examples

``` r
fast_rc("ATCG")
#> [1] "CGAT"
fast_rc(c("ATCG", "AAGG", NA))
#> [1] "CGAT" "CCTT" NA    
fast_rc("AUCG", type = "RNA")
#> [1] "CGAU"
```
