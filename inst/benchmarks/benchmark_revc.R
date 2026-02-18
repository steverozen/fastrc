library(microbenchmark)
library(Biostrings)
library(fastrc)

# 1. Generate 100 sequences of length 30
set.seed(123)
bases <- c("A", "C", "G", "T")
vec_dna <- replicate(
  100,
  paste(sample(bases, 30, replace = TRUE), collapse = "")
)

# 3. Define methods
# Method A: Standard split/paste (wrapped in sapply)
bench_split <- function(v) {
  sapply(v, function(x) {
    chartr("ATGC", "TACG", paste(rev(strsplit(x, NULL)[[1]]), collapse = ""))
  })
}

# Method B: Integer trick (wrapped in sapply)
bench_integer <- function(v) {
  sapply(v, function(x) {
    chartr("ATGC", "TACG", intToUtf8(rev(utf8ToInt(x))))
  })
}

# Method C: Biostrings (Directly vectorized)
bench_bioc <- function(s) {
  Biostrings::reverseComplement(Biostrings::DNAStringSet(s))
}

# Method D: spgs
bench_spgs <- function(s) {
  spgs::reverseComplement(s, content = "dna", case = "as is")
}

bench_seqinr <- function(s) {
  lapply(s, \(x) seqinr::comp(strsplit(x, split = '')))
}

bench_insect <- function(s) {
  insect::rc(s)
}

bench_tktools <- function(s) {
  tktools::revcomp(s)
}

bench_fastrc <- function(s) {
  fastrc::fast_rc(s)
}

# 4. Run Benchmark
results <- microbenchmark(
  split_paste = bench_split(vec_dna),
  integer_trick = bench_integer(vec_dna),
  biostrings = bench_bioc(vec_dna),
  spgs = bench_spgs(vec_dna),
  insect = bench_insect(vec_dna),
  tktools = bench_tktools(vec_dna),
  fastrc = bench_fastrc(vec_dna),
  # siqinr = bench_seqinr(vec_dna),
  times = 100
)

print(results)

vec_dna <- replicate(
  100,
  paste(sample(bases, 10000, replace = TRUE), collapse = "")
)

results <- microbenchmark(
  # split_paste = bench_split(vec_dna),
  # integer_trick = bench_integer(vec_dna),
  biostrings = bench_bioc(vec_dna),
  # spgs = bench_spgs(vec_dna),
  # insect = bench_insect(vec_dna),
  # tktools = bench_tktools(vec_dna),
  # siqinr = bench_seqinr(vec_dna),
  fastrc = bench_fastrc(vec_dna),
  times = 100
)

print(results)


vec_dna <- replicate(
  10,
  paste(sample(bases, 1000000, replace = TRUE), collapse = "")
)

results <- microbenchmark::microbenchmark(
  # split_paste = bench_split(vec_dna),
  #integer_trick = bench_integer(vec_dna),
  biostrings = bench_bioc(vec_dna),
  #spgs = bench_spgs(vec_dna),
  #insect = bench_insect(vec_dna),
  # tktools = bench_tktools(vec_dna),
  # siqinr = bench_seqinr(vec_dna),
  fastrc = bench_fastrc(vec_dna),
  times = 100
)

print(results)
