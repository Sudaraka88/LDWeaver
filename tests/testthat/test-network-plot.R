test_that("create_network adds _ig suffix to intergenic node labels", {
  tophits <- data.frame(
    pos1_ann = c("geneA:anno1", "geneA:anno2", "geneB:anno3", "geneC:anno4"),
    pos2_ann = c("geneB:anno1", "geneC:anno2", "geneC:anno3", "geneD:anno4"),
    MI = c(0.9, 0.8, 0.7, 0.6),
    links = c("nsXig", "igXns", "igXig", "nsXns")
  )

  plot <- LDWeaver::create_network(tophits, min_links_to_include = 1)
  labels <- sort(unique(plot$data$name))

  expect_true("geneA_ig" %in% labels)
  expect_true("geneB_ig" %in% labels)
  expect_true("geneC_ig" %in% labels)
  expect_true("geneD" %in% labels)
  expect_false("geneD_ig" %in% labels)
})
