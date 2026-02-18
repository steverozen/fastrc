#' Fast Reverse Complement of DNA/RNA Sequences
#'
#' Computes the reverse complement of one or more DNA or RNA sequences using a
#' C++ lookup table for O(1) complement mapping. Supports full IUPAC ambiguity
#' codes and preserves case.
#'
#' @param seqs A character vector of DNA or RNA sequences.
#' @param type Either `"DNA"` (A\eqn{\leftrightarrow}T) or `"RNA"`
#'   (A\eqn{\leftrightarrow}U).
#'
#' @return A character vector of reverse-complemented sequences, with `NA`
#'   values preserved.
#'
#' @examples
#' fast_rc("ATCG")
#' fast_rc(c("ATCG", "AAGG", NA))
#' fast_rc("AUCG", type = "RNA")
#'
#' @export
fast_rc <- function(seqs, type = "DNA") {
  if (!is.character(seqs)) {
    stop("`seqs` must be a character vector, not ", class(seqs)[[1]])
  }
  if (!type %in% c("DNA", "RNA")) {
    stop('`type` must be "DNA" or "RNA"')
  }
  fast_rc_impl(seqs, type)
}
