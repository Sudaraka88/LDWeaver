# test_that("Test run on mega sample", {
#   expect_equal(NULL, )
# })

LDWeaver::LDWeaver(dset = "mega",
                   aln_path = system.file("extdata", "sample.aln.gz", package = "LDWeaver"),
                   gbk_path = system.file("extdata", "sample.gbk", package = "LDWeaver"),
                   validate_ref_ann_lengths = F, mega_dset = T)
#
# test_that("Test run on standard sample", {
#   expect_equal(NULL, )
# })

LDWeaver::LDWeaver(dset = "std",
                   aln_path = system.file("extdata", "sample.aln.gz", package = "LDWeaver"),
                   gbk_path = system.file("extdata", "sample.gbk", package = "LDWeaver"),
                   validate_ref_ann_lengths = F)

# path_to_folder = "tests/testthat" # This is for local testing
path_to_folder = getwd()

ths_mega = LDWeaver::read_TopHits(file.path(path_to_folder, "mega/Tophits/sr_tophits.tsv"))
ths_std = LDWeaver::read_TopHits(file.path(path_to_folder, "std/Tophits/sr_tophits.tsv"))

test_that("Short range tophits", {
  expect_equal(all.equal(ths_mega, ths_std), TRUE)
})


thl_mega = LDWeaver::read_TopHits(file.path(path_to_folder, "mega/Tophits/lr_tophits.tsv"))
thl_std = LDWeaver::read_TopHits(file.path(path_to_folder, "std/Tophits/lr_tophits.tsv"))

test_that("Long range tophits", {
  expect_equal(all.equal(thl_mega, thl_std), TRUE)
})

as_mega = LDWeaver::read_AnnotatedLinks(file.path(path_to_folder, "mega/Annotated_links/sr_links_annotated.tsv"))
as_std = LDWeaver::read_AnnotatedLinks(file.path(path_to_folder, "std/Annotated_links/sr_links_annotated.tsv"))

test_that("Short range all links", {
  expect_equal(all.equal(as_mega, as_std), TRUE)
})


al_mega = LDWeaver::read_AnnotatedLinks(file.path(path_to_folder, "mega/Annotated_links/lr_links_annotated.tsv"))
al_std = LDWeaver::read_AnnotatedLinks(file.path(path_to_folder, "std/Annotated_links/lr_links_annotated.tsv"))

test_that("Long range all links", {
  expect_equal(all.equal(al_mega, al_std), TRUE)
})

## Try running on SNP only alignment
pos = as.numeric(readLines(system.file("extdata", "snp_sample.pos", package = "LDWeaver")))
# test_that("Test run on SNP sample", {
#   expect_equal(NULL, )
# })

LDWeaver::LDWeaver(dset = "SNP", aln_has_all_bases = F, pos = pos,
                   aln_path = system.file("extdata", "snp_sample.fa.gz", package = "LDWeaver"),
                   gbk_path = system.file("extdata", "sample.gbk", package = "LDWeaver"),
                   validate_ref_ann_lengths = F, mega_dset = F)
