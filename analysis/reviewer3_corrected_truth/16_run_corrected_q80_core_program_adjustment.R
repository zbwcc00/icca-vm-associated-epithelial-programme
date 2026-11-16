suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(SecAct))

project_root <- "P:"
revision_root <- file.path(project_root, "revision_round1_reviewer3_20261116")
input_root <- file.path(revision_root, "analysis_02_corrected_threshold_pseudobulk", "pseudobulk_inputs")
output_root <- file.path(revision_root, "analysis_07_corrected_q80_core_program_adjustment")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
set.seed(20261116)

definitions <- read.delim(
  file.path(revision_root, "analysis_06_unified_corrected_sensitivities",
            "Table_U0_comparator_programme_gene_definitions.tsv"),
  stringsAsFactors = FALSE
)
programmes <- split(definitions$gene, definitions$program)
programmes <- programmes[c("EMT_Hallmark_anchor_excluded",
                           "Hypoxia_Hallmark_anchor_excluded",
                           "ECM_Reactome_anchor_excluded")]

membership <- read.delim(
  file.path(project_root, "results", "62_secact_receiver_activity",
            "secact_signature_group_membership.tsv"),
  stringsAsFactors = FALSE
)
corrected_external <- read.delim(
  file.path(revision_root, "analysis_04_corrected_secact_integration",
            "Table_R3_13_corrected_q80_source_and_external_constraint.tsv"),
  stringsAsFactors = FALSE
)
corrected_external$external_RNA_compatible <- as.logical(corrected_external$external_RNA_compatible)
corrected_external$nominal_positive_both <- as.logical(corrected_external$nominal_positive_both)
candidate_26 <- corrected_external$signature_group[
  corrected_external$nominal_positive_both & corrected_external$external_RNA_compatible
]
stopifnot(length(unique(candidate_26)) == 26)

fit_one <- function(cohort, adjustment) {
  d <- file.path(input_root, paste0(cohort, "_q80"))
  counts <- as(readMM(file.path(d, "pseudobulk_counts_genes_by_samples.mtx")), "dgCMatrix")
  genes <- read.delim(file.path(d, "genes.tsv"), stringsAsFactors = FALSE)
  samples <- read.delim(file.path(d, "pseudobulk_samples.tsv"), stringsAsFactors = FALSE)
  rownames(counts) <- make.unique(genes$gene_symbol)
  colnames(counts) <- samples$pseudobulk_id
  samples$replicate_id <- factor(samples$library_id)
  samples$score_group <- factor(samples$score_group, levels = c("low", "high"))

  dge_for_score <- calcNormFactors(DGEList(counts = counts))
  logcpm <- cpm(dge_for_score, log = TRUE, prior.count = 0.5)
  if (adjustment != "unadjusted") {
    present <- intersect(programmes[[adjustment]], rownames(logcpm))
    if (length(present) < 5) stop("Too few genes for ", adjustment)
    gene_z <- t(scale(t(logcpm[present, , drop = FALSE])))
    gene_z[!is.finite(gene_z)] <- 0
    raw <- colMeans(gene_z)
    centered <- raw - ave(raw, samples$replicate_id, FUN = mean)
    samples$within_pair_programme <- centered / sd(centered)
    design <- model.matrix(~ replicate_id + score_group + within_pair_programme, data = samples)
  } else {
    design <- model.matrix(~ replicate_id + score_group, data = samples)
  }
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient design")
  coefficient <- which(colnames(design) == "score_grouphigh")
  y <- DGEList(counts = counts)
  keep <- filterByExpr(y, design = design, min.count = 10)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = coefficient)
  differential <- topTags(test, n = Inf, sort.by = "none")$table
  profile <- matrix(differential$logFC, ncol = 1,
                    dimnames = list(rownames(differential), "high_minus_remainder"))
  result <- SecAct.activity.inference(
    inputProfile = profile, is.differential = TRUE, sigMatrix = "SecAct",
    is.filter.sig = FALSE, is.group.sig = TRUE, lambda = 5e5,
    nrand = 1000, ncores = 1L, backend = "cpu-pure", rng_method = "mt19937"
  )
  protein <- data.frame(
    secreted_protein = rownames(result$zscore), beta = result$beta[, 1],
    zscore = result$zscore[, 1], pvalue = result$pvalue[, 1],
    stringsAsFactors = FALSE
  )
  protein <- merge(protein, membership, by = "secreted_protein", all.x = TRUE)
  group <- aggregate(protein[, c("beta", "zscore", "pvalue")],
                     by = protein[, c("signature_group", "group_members", "group_size")],
                     FUN = function(x) x[1])
  group$cohort <- cohort
  group$adjustment <- adjustment
  group$design_condition_number <- kappa(design)
  group$residual_df_design <- nrow(design) - ncol(design)
  group
}

cohorts <- c("GSE138709", "GSE125449")
adjustments <- c("unadjusted", names(programmes))
fits <- list()
for (cohort in cohorts) {
  for (adjustment in adjustments) {
    message(cohort, " / ", adjustment)
    fits[[paste(cohort, adjustment, sep = "__")]] <- fit_one(cohort, adjustment)
  }
}
all_results <- do.call(rbind, fits)
rownames(all_results) <- NULL
write.table(all_results,
            gzfile(file.path(output_root, "corrected_q80_core_program_adjustment_all_groups.tsv.gz"), "wt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

rows <- list()
for (cohort in cohorts) {
  reference <- all_results[all_results$cohort == cohort & all_results$adjustment == "unadjusted", ]
  for (adjustment in names(programmes)) {
    adjusted <- all_results[all_results$cohort == cohort & all_results$adjustment == adjustment, ]
    merged <- merge(reference[, c("signature_group", "zscore", "beta")],
                    adjusted[, c("signature_group", "zscore", "beta", "pvalue",
                                 "design_condition_number", "residual_df_design")],
                    by = "signature_group", suffixes = c("_unadjusted", "_adjusted"))
    for (set_name in c("all_657", "corrected_q80_external_26")) {
      keep <- if (set_name == "all_657") rep(TRUE, nrow(merged)) else merged$signature_group %in% candidate_26
      x <- merged[keep, ]
      rows[[length(rows) + 1]] <- data.frame(
        cohort = cohort, adjustment = adjustment, candidate_set = set_name,
        groups = nrow(x),
        spearman_zscore_unadjusted_vs_adjusted = cor(x$zscore_unadjusted, x$zscore_adjusted, method = "spearman"),
        same_beta_direction_fraction = mean(sign(x$beta_unadjusted) == sign(x$beta_adjusted)),
        positive_adjusted = sum(x$beta_adjusted > 0),
        nominal_positive_adjusted = sum(x$beta_adjusted > 0 & x$pvalue < .05),
        design_condition_number = unique(x$design_condition_number),
        residual_df_design = unique(x$residual_df_design),
        stringsAsFactors = FALSE
      )
    }
  }
}
summary <- do.call(rbind, rows)
write.table(summary, file.path(output_root, "Table_U11_corrected_q80_core_program_adjustment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
print(summary)
