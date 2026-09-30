# Regression tests for bugs found in the package audit

e12 <- rbind(c("A","G"), c("G","A"), c("T","C"), c("C","T"), c("A","T"), c("T","A"),
             c("G","T"), c("T","G"), c("G","C"), c("C","G"), c("A","C"), c("C","A"))

# Genomic A site edited 75% in treat and 30% in control: G is the pooled majority
make_heavy_site <- function() {
  mk <- function(a, g) data.frame(chr = "chr1", pos = 100L, ref = "A",
                                  A = a, T = 0, C = 0, G = g, site_id = "chr1_100")
  list(t1 = mk(25, 75), t2 = mk(25, 75), t3 = mk(25, 75),
       c1 = mk(70, 30), c2 = mk(70, 30))
}
heavy_sig <- data.frame(
  groupID = "chr1_100", featureID = c("EG", "EA"), padj = 1e-5,
  pvalue = c(1e-6, 2e-6), log2fold_treat_control = c(2, -1),
  control = 0, treat = 0, stringsAsFactors = FALSE
)
heavy_design <- c(t1 = "treat", t2 = "treat", t3 = "treat",
                  c1 = "control", c2 = "control")

test_that("heavily edited site keeps genome ref instead of pooled majority", {
  res <- compute_site_meta("chr1_100", heavy_sig, make_heavy_site(), e12,
                           heavy_design, c("t1", "t2", "t3"), c("c1", "c2"),
                           symmetric = FALSE, genome_ref = "A")
  expect_equal(res$ref, "A")
  expect_equal(res$targ, "G")
  expect_false(res$snp)
  expect_equal(res$prop, 0.75)
  expect_equal(res$fold_change, 2)
})

test_that("SNP vs genome: ref is the base carried by the controls, flagged as snp", {
  # Genome says G, but all lines carry A; low A>G editing in treat
  mk <- function(a, g) data.frame(A = a, T = 0, C = 0, G = g, site_id = "s_1")
  dl <- list(t1 = mk(90, 10), t2 = mk(92, 8), c1 = mk(99, 1), c2 = mk(100, 0))
  sig <- data.frame(groupID = "s_1", featureID = c("EG", "EA"), padj = 1e-4,
                    pvalue = c(1e-5, 2e-5), log2fold_treat_control = c(3, -0.1),
                    control = 0, treat = 0)
  res <- compute_site_meta("s_1", sig, dl, e12, NULL, c("t1", "t2"), c("c1", "c2"),
                           symmetric = FALSE, genome_ref = "G")
  expect_equal(res$ref, "A")
  expect_equal(res$targ, "G")
  expect_true(res$snp)
  expect_equal(res$prop, 18 / 200)
})

test_that("infer_site_ref resolves control and treat majorities that disagree", {
  m <- rbind(t1 = c(A = 20, T = 0, C = 0, G = 80), c1 = c(A = 90, T = 0, C = 0, G = 10))
  # Control majority wins when both agree with / without genome
  expect_equal(infer_site_ref(m, "t1", "c1", NA, e12), "A")
  # Control edited above 50% (e.g. FvF): genome base decides
  m2 <- m[c("c1", "t1"), ]; rownames(m2) <- c("t1", "c1")
  expect_equal(infer_site_ref(m2, "t1", "c1", "A", e12), "A")
  # ... and without genome base, the valid edit direction decides
  expect_equal(infer_site_ref(m2, "t1", "c1", NA, rbind(c("A", "G"))), "A")
  # No control coverage: fall back to treat majority, then genome
  m3 <- rbind(t1 = c(A = 9, T = 0, C = 0, G = 1), c1 = c(A = 0, T = 0, C = 0, G = 0))
  expect_equal(infer_site_ref(m3, "t1", "c1", NA, e12), "A")
  expect_equal(infer_site_ref(m3 * 0, "t1", "c1", "C", e12), "C")
})

test_that("get_hits takes the genome ref from the data when ref_gr is NULL", {
  hits <- get_hits(heavy_sig, fdr = 0.1, n_cores = 1, data_list = make_heavy_site(),
                   edits_of_interest = e12, design_vector = heavy_design)
  expect_length(hits, 1)
  expect_equal(hits$ref, "A")
  expect_equal(hits$targ, "G")
  expect_equal(hits$base, "A")
  expect_false(hits$snp)
})

test_that("get_hits takes the genome ref from ref_gr when given", {
  ref_gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, width = 1), ref = "A")
  dl <- lapply(make_heavy_site(), function(x) { x$ref <- NULL; x })
  hits <- get_hits(heavy_sig, fdr = 0.1, n_cores = 1, data_list = dl,
                   edits_of_interest = e12, design_vector = heavy_design,
                   ref_gr = ref_gr)
  expect_equal(hits$ref, "A")
  expect_equal(hits$targ, "G")
})

test_that("site with no valid target does not crash get_hits", {
  # Dominant/genome C, significant EG and EC, only A>G and T>C of interest
  dl <- list(
    t1 = data.frame(ref = "C", A = 0, T = 0, C = 90, G = 10, site_id = "x_1"),
    c1 = data.frame(ref = "C", A = 0, T = 0, C = 99, G = 1,  site_id = "x_1")
  )
  sig <- data.frame(groupID = "x_1", featureID = c("EG", "EC"), padj = 0.01,
                    pvalue = c(0.001, 0.002), log2fold_treat_control = 1,
                    control = 0, treat = 0)
  dv <- c(t1 = "treat", c1 = "control")
  hits <- get_hits(sig, fdr = 0.1, n_cores = 1, data_list = dl, design_vector = dv)
  # EC is the ref itself; EG is kept as a (non-ADAR) C>G hit ...
  expect_equal(hits$targ, "G")
  # ... which the ADAR edit filter then removes
  expect_length(suppressMessages(filter_edits(list(e = hits)))$e, 0)

  # Only the ref base itself significant: no edit left, site dropped
  hits_ref_only <- get_hits(sig[2, ], fdr = 0.1, n_cores = 1, data_list = dl,
                            design_vector = dv)
  expect_length(hits_ref_only, 0)
})

test_that("symmetric mode tolerates NA fold change", {
  sig <- heavy_sig[1, ]
  sig$log2fold_treat_control <- NA_real_
  expect_no_error(
    compute_site_meta("chr1_100", sig, make_heavy_site(), e12, heavy_design,
                      c("t1", "t2", "t3"), c("c1", "c2"), symmetric = TRUE,
                      genome_ref = "A")
  )
})

test_that("flanked site spanning CDS and UTR is annotated without error", {
  gtf <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(100, 100, 100, 161), c(200, 200, 160, 200)),
    strand = "+", type = c("exon", "CDS", "CDS", "3UTR"),
    gene_id = "G", transcript_id = "T1"
  )
  site <- GenomicRanges::GRanges("chr1", IRanges::IRanges(250, width = 1), strand = "+")
  res <- annotate_with_genes(site, gtf, c(G = "g"), c(T1 = 1))
  expect_equal(unname(res$gene), "G")
  expect_true(res$out_of_range)
  expect_equal(unname(res$transcript_types), "CDS")
})

test_that("site in a stop codon is annotated without error", {
  gtf <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(100, 198), c(200, 200)), strand = "+",
    type = c("exon", "stop_codon"), gene_id = "G", transcript_id = "T1"
  )
  site <- GenomicRanges::GRanges("chr1", IRanges::IRanges(199, width = 1))
  res <- annotate_with_genes(site, gtf, c(G = "g"), c(T1 = 1))
  expect_equal(unname(res$transcript_types), "stop_codon")
})

test_that("addStrandForHyperTRIBE sets strand from the ADAR edit type", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1:3, width = 1))
  gr$ref  <- c("A", "T", "C")
  gr$targ <- c("G", "C", "T")
  expect_equal(as.character(BiocGenerics::strand(addStrandForHyperTRIBE(gr))),
               c("+", "-", "*"))
})

test_that("annotate_hits does not assign T>C sites to + strand genes", {
  gtf <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(100, 1e5), c(200, 1e5 + 100)),
    strand = c("+", "-"), type = "exon",
    gene_id = c("PLUS", "MINUS"), transcript_id = c("TP", "TM")
  )
  sites <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(150, 160, 1e5 + 50), width = 1)
  )
  sites$ref  <- c("A", "T", "T")
  sites$targ <- c("G", "C", "C")
  design <- tibble::tibble(experiment = "e1", control = list("c1"))
  txi <- list(abundance = matrix(1, nrow = 2, ncol = 1,
                                 dimnames = list(c("TP", "TM"), "c1")))

  res <- annotate_hits(list(e1 = sites), design, txi, gtf,
                       c(PLUS = "plus", MINUS = "minus"))$e1

  expect_equal(unname(res$gene[1]), "PLUS")
  expect_true(is.na(res$gene[2]))        # T>C inside a + gene: no + gene match
  expect_equal(unname(res$gene[3]), "MINUS")
  expect_equal(as.character(BiocGenerics::strand(res)), c("+", "-", "-"))
})

test_that("filter_edits drops NaN proportions instead of erroring", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1:3, width = 1))
  gr$ref  <- "A"
  gr$targ <- "G"
  gr$prop <- c(0.2, NaN, 0.5)
  gr$fold_change <- c(2, 2, NA)
  res <- suppressMessages(filter_edits(list(e1 = gr), fc_thresh = 1))$e1
  expect_equal(BiocGenerics::start(res), 1L)
})

test_that("stranded extract_count_data complements ref on reverse rows", {
  # one sample: A T C G a t c g
  dat <- data.frame(V1 = "chr1", V2 = 10L, V3 = "a",
                    V4 = 5L, V5 = 0L, V6 = 0L, V7 = 1L,
                    V8 = 0L, V9 = 7L, V10 = 2L, V11 = 0L)
  res <- suppressMessages(extract_count_data(dat, "s1", stranded = TRUE))$s1
  fwd <- res[res$strand == "+", ]
  rev <- res[res$strand == "-", ]
  expect_equal(fwd$ref, "A")
  expect_equal(rev$ref, "T")
  # reverse A/T/C/G are the complemented lowercase t/a/g/c counts
  expect_equal(unlist(rev[, c("A", "T", "C", "G")]), c(A = 7, T = 0, C = 0, G = 2))
})

test_that("generate_count_files removes stale count files", {
  tmpdir <- withr::local_tempdir()
  file.create(file.path(tmpdir, "counts_Treat_OLD.txt"))
  dl <- list(
    T1 = data.frame(chr = "c", pos = 1L, ref = "A", A = 5L, T = 0L, C = 0L, G = 1L, site_id = "c_1"),
    C1 = data.frame(chr = "c", pos = 1L, ref = "A", A = 6L, T = 0L, C = 0L, G = 0L, site_id = "c_1")
  )
  suppressMessages(generate_count_files(dl, c("treat", "control"), tmpdir))
  expect_setequal(list.files(tmpdir, pattern = "^counts_"),
                  c("counts_Treat_T1.txt", "counts_Control_C1.txt"))
})
