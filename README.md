# fastrc

Fast reverse complement of DNA and RNA sequences in R, implemented in C++ via Rcpp.

`fastrc` uses a static lookup table for O(1) per-base complement mapping with full IUPAC ambiguity code support. It is **especially useful for reverse complementing many short sequences** (e.g. primers, probes, k-mers, short reads), where per-call overhead dominates and `fastrc` is over 100x faster than alternatives.

## Installation

```r
# From GitHub
devtools::install_github("steverozen/fastrc")
```

## Usage

```r
library(fastrc)

fast_rc("ATCG")
#> [1] "CGAT"

fast_rc(c("ATCG", "AAGG", NA))
#> [1] "CGAT" "CCTT" NA

fast_rc("AUCG", type = "RNA")
#> [1] "CGAU"
```

## Features

- DNA (A&#8596;T) and RNA (A&#8596;U) modes
- Full IUPAC ambiguity code support (M&#8596;K, R&#8596;Y, S&#8596;S, W&#8596;W, V&#8596;B, H&#8596;D, N&#8596;N)
- Case preservation
- NA handling
- Vectorized over character vectors

## Benchmarks

Benchmarks were run on a 12th Gen Intel i7-1270P. The package was compiled
with `-O3 -flto` (the default for installed packages).

### 100 sequences x 30 bp

This is the scenario where `fastrc` shines brightest -- many short sequences
where per-call overhead matters most.

| Method | Median | vs fastrc |
|---|---:|---:|
| **fastrc** | **13 &#181;s** | **1x** |
| spgs | 621 &#181;s | 48x slower |
| integer trick | 675 &#181;s | 52x slower |
| insect | 1,032 &#181;s | 79x slower |
| Biostrings | 1,579 &#181;s | 121x slower |
| tktools | 3,351 &#181;s | 258x slower |
| split/paste | 4,375 &#181;s | 337x slower |

### 100 sequences x 10,000 bp

| Method | Median | vs fastrc |
|---|---:|---:|
| **fastrc** | **1.4 ms** | **1x** |
| Biostrings | 3.3 ms | 2.3x slower |

### 10 sequences x 1,000,000 bp

| Method | Median | vs fastrc |
|---|---:|---:|
| **fastrc** | **14.2 ms** | **1x** |
| Biostrings | 24.0 ms | 1.7x slower |

## License

GPL (>= 3)
