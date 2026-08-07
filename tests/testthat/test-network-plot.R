test_that("create_network adds _ig suffix to intergenic node labels (include_intergenic = TRUE)", {
  tophits <- data.frame(
    pos1_ann = c("geneA:anno1", "geneA:anno2", "geneB:anno3", "geneC:anno4"),
    pos2_ann = c("geneB:anno1", "geneC:anno2", "geneC:anno3", "geneD:anno4"),
    MI = c(0.9, 0.8, 0.7, 0.6),
    links = c("nsXig", "igXns", "igXig", "nsXns")
  )

  plot <- LDWeaver::create_network(tophits, min_links_to_include = 1, include_intergenic = TRUE)
  labels <- sort(unique(plot$data$name))

  # intergenic sides get _ig suffix
  expect_true("geneA_ig" %in% labels)
  expect_true("geneB_ig" %in% labels)
  expect_true("geneC_ig" %in% labels)
  # genic-only node (nsXns) appears without _ig
  expect_true("geneD" %in% labels)
  expect_false("geneD_ig" %in% labels)
  # geneA also appears genic-only (pos1 of igXns is geneA with ig, but here nsXig has pos1 as geneA genic)
  expect_true("geneA" %in% labels)
})

test_that("create_network with include_intergenic = FALSE filters out ig links", {
  tophits <- data.frame(
    pos1_ann = c("geneA:anno1", "geneA:anno2", "geneB:anno3", "geneC:anno4"),
    pos2_ann = c("geneB:anno1", "geneC:anno2", "geneC:anno3", "geneD:anno4"),
    MI = c(0.9, 0.8, 0.7, 0.6),
    links = c("nsXig", "igXns", "igXig", "nsXns")
  )

  plot <- LDWeaver::create_network(tophits, min_links_to_include = 1, include_intergenic = FALSE)
  labels <- sort(unique(plot$data$name))

  # only the nsXns link should remain; no _ig labels
  expect_false(any(grepl("_ig$", labels)))
  expect_true("geneC" %in% labels)
  expect_true("geneD" %in% labels)
})

test_that("create_network returns NULL with a message when no genic links exist", {
  tophits <- data.frame(
    pos1_ann = c("geneA:anno1", "geneB:anno2"),
    pos2_ann = c("geneC:anno1", "geneD:anno2"),
    MI = c(0.9, 0.8),
    links = c("nsXig", "igXns")
  )

  expect_message(
    result <- LDWeaver::create_network(tophits, min_links_to_include = 1, include_intergenic = FALSE),
    "No genic-only links available"
  )
  expect_null(result)
})
