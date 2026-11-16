suppressPackageStartupMessages(library(Matrix))
loadNamespace("edgeR", lib.loc = "P:/.R_library")
loadNamespace("SecAct", lib.loc = "P:/.R_library")

root <- "P:/revision_round1_reviewer3_20261116"
input_root <- file.path(root, "analysis_02_corrected_threshold_pseudobulk", "pseudobulk_inputs")
output_root <- file.path(root, "analysis_03_corrected_edger_secact")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

run_one <- function(input_dir) {
  label <- basename(input_dir)
  counts <- as(readMM(file.path(input_dir, "pseudobulk_counts_genes_by_samples.mtx")), "dgCMatrix")
  genes <- read.delim(file.path(input_dir, "genes.tsv"), check.names = FALSE)
  metadata <- read.delim(file.path(input_dir, "pseudobulk_samples.tsv"), check.names = FALSE)
  rownames(counts) <- make.unique(genes$gene_symbol)
  colnames(counts) <- metadata$pseudobulk_id
  metadata$score_group <- factor(metadata$score_group, levels = c("low", "high"))
  metadata$library_id <- factor(metadata$library_id)

  y <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(y, group = metadata$score_group, min.count = 10)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  design <- model.matrix(~ library_id + score_group, data = metadata)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  test <- edgeR::glmQLFTest(fit, coef = "score_grouphigh")
  differential <- edgeR::topTags(test, n = Inf, sort.by = "PValue")$table
  differential$gene_symbol <- rownames(differential)
  differential <- differential[, c("gene_symbol", setdiff(colnames(differential), "gene_symbol"))]
  write.table(differential, file.path(output_root, paste0(label, "_high_vs_low_edgeR.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  profile <- matrix(differential$logFC, ncol = 1,
                    dimnames = list(differential$gene_symbol, "high_minus_low"))
  set.seed(20261116)
  result <- SecAct::SecAct.activity.inference(
    inputProfile = profile,
    is.differential = TRUE,
    sigMatrix = "P:/.R_library/SecAct/extdata/SecAct.tsv.gz",
    is.filter.sig = FALSE,
    is.group.sig = TRUE,
    lambda = 5e5,
    nrand = 1000,
    ncores = 1L,
    backend = "cpu-pure",
    rng_method = "mt19937"
  )
  for (metric in c("beta", "se", "zscore", "pvalue")) {
    values <- result[[metric]]
    output <- data.frame(secreted_protein = rownames(values), value = values[, 1],
                         row.names = NULL, check.names = FALSE)
    colnames(output)[2] <- metric
    write.table(output, file.path(output_root, paste0(label, "_SecAct_", metric, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  data.frame(
    analysis = label,
    libraries = length(unique(metadata$library_id)),
    patients = length(unique(metadata$patient_id)),
    pseudobulks = nrow(metadata),
    genes_tested_edgeR = nrow(differential),
    SecAct_expanded_members = nrow(result$zscore),
    stringsAsFactors = FALSE
  )
}

input_dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
summaries <- do.call(rbind, lapply(input_dirs, run_one))
write.table(summaries, file.path(output_root, "corrected_threshold_edgeR_SecAct_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
session_lines <- c(
  capture.output(R.version),
  "edgeR=4.4.2",
  "SecAct=1.1.0",
  "Package versions were read from installed DESCRIPTION files because sessionInfo() could not traverse the non-ASCII project library path under this Windows locale."
)
writeLines(session_lines, file.path(output_root, "sessionInfo.txt"))
writeLines(c(
  "edgeR design=~ library_id + score_group",
  "score_group=corrected six-gene high-score versus low-score operational split",
  "thresholds=within-cohort 80th and 90th percentiles; endothelial exclusion fixed at 90th percentile",
  "SecAct nrand=1000; seed=20261116; backend=cpu-pure; lambda=5e5",
  "claim boundary=threshold sensitivity of receiver-side transcript compatibility only; not secretion, receptor engagement, communication or causality"
), file.path(output_root, "analysis_manifest.txt"))
