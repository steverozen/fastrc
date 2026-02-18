# fastrc

Fast reverse complement of DNA and RNA sequences in R, implemented in C++ via Rcpp.

`fastrc` uses a static lookup table for O(1) per-base complement mapping with full IUPAC ambiguity code support. It is **especially useful for reverse complementing many short sequences** (e.g. primers, probes, k-mers, short reads), where per-call overhead dominates and `fastrc` is nearly 100x faster than the implementation
in `Biostrings`.

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

Benchmarks were run on a 12th Gen Intel i7-1270P using R's default
compilation flags (`-O2`). See `inst/benchmarks/benchmark_revc.R` to
reproduce.

### 100 sequences x 30 bp

This is where you may really want to use `fastrc`: many short sequences
where per-call overhead matters most.

| Method | Median | vs fastrc |
|---|---:|---:|
| **fastrc** | **17 &#181;s** | **1x** |
| spgs | 619 &#181;s | 36x slower |
| insect | 1,031 &#181;s | 61x slower |
| Biostrings | 1,613 &#181;s | 95x slower |
| tktools | 3,277 &#181;s | 193x slower |

### 100 sequences x 10,000 bp

| Method | Median | vs fastrc |
|---|---:|---:|
| **fastrc** | **1.4 ms** | **1x** |
| Biostrings | 3.4 ms | 2.4x slower |

### 10 sequences x 1,000,000 bp

| Method | Median | vs fastrc |
|---|---:|---:|
| **fastrc** | **14.2 ms** | **1x** |
| Biostrings | 23.8 ms | 1.7x slower |

### Faster local builds

You can get additional speed by adding optimization flags to your
`~/.R/Makevars` file:

```
CXXFLAGS += -O3 -march=native -flto
```

Then reinstall the package. This typically yields a 10-20% improvement
on longer sequences.

## License

GPL (>= 3)
