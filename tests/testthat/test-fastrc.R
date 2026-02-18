test_that("basic DNA reverse complement works", {
  expect_equal(fast_rc("ATCG"), "CGAT")
  expect_equal(fast_rc("AAAA"), "TTTT")
  expect_equal(fast_rc("CCCC"), "GGGG")
  expect_equal(fast_rc("TTTT"), "AAAA")
  expect_equal(fast_rc("GGGG"), "CCCC")
})

test_that("single base complements are correct", {
  expect_equal(fast_rc("A"), "T")
  expect_equal(fast_rc("T"), "A")
  expect_equal(fast_rc("C"), "G")
  expect_equal(fast_rc("G"), "C")
})

test_that("reverse complement reverses order", {
  expect_equal(fast_rc("ACGT"), "ACGT")
  expect_equal(fast_rc("AACG"), "CGTT")
})

test_that("case is preserved", {
  expect_equal(fast_rc("atcg"), "cgat")
  expect_equal(fast_rc("AtCg"), "cGaT")
})

test_that("RNA mode uses U instead of T", {
  expect_equal(fast_rc("AUCG", type = "RNA"), "CGAU")
  expect_equal(fast_rc("A", type = "RNA"), "U")
  expect_equal(fast_rc("U", type = "RNA"), "A")
  expect_equal(fast_rc("a", type = "RNA"), "u")
})

test_that("IUPAC ambiguity codes are complemented", {
  expect_equal(fast_rc("M"), "K")
  expect_equal(fast_rc("K"), "M")
  expect_equal(fast_rc("R"), "Y")
  expect_equal(fast_rc("Y"), "R")
  expect_equal(fast_rc("S"), "S")
  expect_equal(fast_rc("W"), "W")
  expect_equal(fast_rc("V"), "B")
  expect_equal(fast_rc("B"), "V")
  expect_equal(fast_rc("H"), "D")
  expect_equal(fast_rc("D"), "H")
  expect_equal(fast_rc("N"), "N")
})

test_that("lowercase IUPAC codes are complemented", {
  expect_equal(fast_rc("m"), "k")
  expect_equal(fast_rc("k"), "m")
  expect_equal(fast_rc("r"), "y")
  expect_equal(fast_rc("y"), "r")
  expect_equal(fast_rc("v"), "b")
  expect_equal(fast_rc("b"), "v")
  expect_equal(fast_rc("h"), "d")
  expect_equal(fast_rc("d"), "h")
  expect_equal(fast_rc("n"), "n")
})

test_that("NA values are preserved", {
  result <- fast_rc(c("ATCG", NA, "GGCC"))
  expect_equal(result[1], "CGAT")
  expect_true(is.na(result[2]))
  expect_equal(result[3], "GGCC")
})

test_that("vectorized input works", {
  input <- c("ATCG", "AAGG", "TTCC")
  expected <- c("CGAT", "CCTT", "GGAA")
  expect_equal(fast_rc(input), expected)
})

test_that("empty string returns empty string", {
  expect_equal(fast_rc(""), "")
})

test_that("length-30 sequences are reverse complemented correctly", {
  input <- c("ATCGATCGATCGATCGATCGATCGATCGAT",
             "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGC")
  expected <- c("ATCGATCGATCGATCGATCGATCGATCGAT",
                "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGC")
  expect_equal(fast_rc(input), expected)
  expect_equal(nchar(input[1]), 30)
  expect_equal(nchar(input[2]), 30)
})

test_that("non-character input raises an error", {
  expect_error(fast_rc(42), "must be a character vector")
  expect_error(fast_rc(list("ATCG")), "must be a character vector")
  expect_error(fast_rc(factor("ATCG")), "must be a character vector")
})

test_that("invalid type raises an error", {
  expect_error(fast_rc("ATCG", type = "PROTEIN"), 'must be "DNA" or "RNA"')
})

test_that("switching between DNA and RNA modes works", {
  expect_equal(fast_rc("A", type = "DNA"), "T")
  expect_equal(fast_rc("A", type = "RNA"), "U")
  expect_equal(fast_rc("A", type = "DNA"), "T")
})
