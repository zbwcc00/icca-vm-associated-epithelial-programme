project_root <- "P:"
revision_root <- file.path(project_root, "revision_round1_reviewer3_20261116")
output_root <- file.path(revision_root, "analysis_09_corrected_secact_redundancy")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(SecAct))

signature <- SecAct:::load_sig_matrix("SecAct", 5e5)$X
grouped <- as.matrix(SecAct:::cluster_signatures(signature, is.group.cor = 0.9))
mapping <- unique(read.delim(
  file.path(project_root, "results", "62_secact_receiver_activity", "secact_signature_group_membership.tsv"),
  stringsAsFactors = FALSE
)[, c("signature_group", "group_members", "group_size")])
mapping <- mapping[order(mapping$signature_group), ]
if (ncol(grouped) != nrow(mapping)) stop("Grouped signature count mismatch")
colnames(grouped) <- mapping$signature_group
correlation <- cor(grouped, method = "pearson")

cross <- read.delim(
  file.path(revision_root, "analysis_04_corrected_secact_integration", "corrected_cross_cohort_SecAct_q80.tsv"),
  stringsAsFactors = FALSE
)
cross$nominal_positive_both <- as.logical(cross$nominal_positive_both)
cross$global_FDR_positive_both <- as.logical(cross$global_FDR_positive_both)
external <- read.delim(
  file.path(revision_root, "analysis_04_corrected_secact_integration", "Table_R3_13_corrected_q80_source_and_external_constraint.tsv"),
  stringsAsFactors = FALSE
)
external$has_source_restricted_member <- as.logical(external$has_source_restricted_member)
external$external_RNA_compatible <- as.logical(external$external_RNA_compatible)
external$nominal_positive_both <- as.logical(external$nominal_positive_both)

effective_rank <- function(m) {
  values <- pmax(eigen(m, symmetric = TRUE, only.values = TRUE)$values, 0)
  p <- values/sum(values); p <- p[p > 0]
  c(participation_ratio = sum(values)^2/sum(values^2),
    shannon_effective_rank = exp(-sum(p*log(p))))
}

summarize_set <- function(name, ids) {
  ids <- intersect(ids, colnames(correlation))
  m <- correlation[ids, ids, drop = FALSE]
  pairs <- m[lower.tri(m)]
  rank <- effective_rank(m)
  data.frame(
    set = name, groups = length(ids), pair_count = length(pairs),
    median_correlation = median(pairs), median_absolute_correlation = median(abs(pairs)),
    fraction_abs_r_ge_0_5 = mean(abs(pairs) >= .5),
    fraction_abs_r_ge_0_7 = mean(abs(pairs) >= .7),
    fraction_r_ge_0_8 = mean(pairs >= .8), fraction_r_ge_0_9 = mean(pairs >= .9),
    participation_ratio = unname(rank["participation_ratio"]),
    shannon_effective_rank = unname(rank["shannon_effective_rank"]),
    stringsAsFactors = FALSE
  )
}

sets <- list(
  all_657_groups = mapping$signature_group,
  corrected_q80_nominal_268 = cross$signature_group[cross$nominal_positive_both],
  corrected_q80_global_FDR_265 = cross$signature_group[cross$global_FDR_positive_both],
  corrected_q80_source_restricted_58 = external$signature_group[
    external$nominal_positive_both & external$has_source_restricted_member
  ],
  corrected_q80_external_26 = external$signature_group[
    external$nominal_positive_both & external$external_RNA_compatible
  ]
)
stopifnot(length(unique(sets$corrected_q80_source_restricted_58)) == 58)
stopifnot(length(unique(sets$corrected_q80_external_26)) == 26)
summary <- do.call(rbind, Map(summarize_set, names(sets), sets))
write.table(summary, file.path(output_root, "Table_U13_corrected_SecAct_signature_redundancy.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  "Corrected q80 SecAct target-signature redundancy audit.",
  "The 657 correlation-collapsed groups are correlated hypotheses, not independent discoveries.",
  "The 26-group set is defined by the corrected q80 external source-RNA compatibility funnel."
), file.path(output_root, "README.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
print(summary)
