test_that("read_gene_list accepts a character vector", {
  genes <- read_gene_list(c("TP53", "KRAS", "TP53", " BRAF "))
  expect_equal(genes, c("TP53", "KRAS", "BRAF"))
})

test_that("read_gene_list treats a single non-file string as one gene", {
  expect_equal(read_gene_list("TP53"), "TP53")
})

test_that("read_gene_list reads the legacy two-column CSV format", {
  path <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(data.frame(x = c("TP53", "KRAS", "APC")), path)
  expect_equal(read_gene_list(path), c("TP53", "KRAS", "APC"))
})

test_that("read_gene_list reads a newline-delimited txt file", {
  path <- withr::local_tempfile(fileext = ".txt")
  writeLines(c("TP53", "KRAS", "", "APC"), path)
  expect_equal(read_gene_list(path), c("TP53", "KRAS", "APC"))
})
