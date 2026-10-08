#####################################################################
# QC check: Effect of batch covariate on GWAS calibration (Factor 1)
# Compares lambda_GC and per-SNP p-value concordance with vs without
# batch as a covariate, and checks whether large per-SNP shifts occur
# among significant or non-significant SNPs.
#####################################################################

# interactive -t 5:00:00  -A uppmax2025-2-420
# module load R-bundle-CRAN
# R

library(dplyr)
library(ggplot2)

#####################################################################
# 1. Load GWAS summary stats
#    Expected columns (adjust names as needed): SNP, CHR, BP, P
#    e.g. GCTA --mlma output: Chr SNP bp A1 A2 Freq b se p
#####################################################################

setwd("/proj/uppmax2025-2-49/snic2022-6-229/private/Darwin_Affy/NEW_DA_IMP_MERGED/2024-10-11/")
gwas_noBatch <- read.table("CCDF1_QC6_sex_age.mlma",       header = TRUE, stringsAsFactors = FALSE)
gwas_Batch   <- read.table("CCDF1_QC6_sex_age_batch.mlma", header = TRUE, stringsAsFactors = FALSE)

# Standardize column names if needed (edit mapping to match your files)
std_cols <- function(df) {
  names(df) <- tolower(names(df))
  df <- df %>% rename(chr = any_of(c("chr", "chromosome")),
                       bp  = any_of(c("bp", "pos", "position")),
                       p   = any_of(c("p", "pval", "p_value")))
  df
}
gwas_noBatch <- std_cols(gwas_noBatch)
gwas_Batch   <- std_cols(gwas_Batch)

#####################################################################
# 2. Lambda_GC calculation
#####################################################################

lambda_gc <- function(pvals) {
  pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals <= 1]
  chisq <- qchisq(1 - pvals, df = 1)
  median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
}

lambda_noBatch <- lambda_gc(gwas_noBatch$p)
lambda_Batch   <- lambda_gc(gwas_Batch$p)

cat("==================== LAMBDA_GC COMPARISON ====================\n")
cat(sprintf("Model WITHOUT batch (age + sex only): lambda_GC = %.4f\n", lambda_noBatch))
cat(sprintf("Model WITH batch (age + sex + batch): lambda_GC = %.4f\n", lambda_Batch))
cat(sprintf("Difference (no_batch - with_batch):   %.4f\n", lambda_noBatch - lambda_Batch))
cat("================================================================\n\n")

#####################################################################
# 3. Per-SNP p-value correlation and magnitude of change
#    Compares -log10(p) between runs for all SNPs present in both
#####################################################################

pval_compare <- inner_join(
  gwas_noBatch %>% select(chr, bp, p) %>% rename(p_noBatch = p),
  gwas_Batch   %>% select(chr, bp, p) %>% rename(p_Batch = p),
  by = c("chr", "bp")) %>%
  mutate(
    log10p_noBatch = -log10(p_noBatch),
    log10p_Batch   = -log10(p_Batch),
    log10p_diff    = log10p_Batch - log10p_noBatch,
    abs_diff       = abs(log10p_diff),
    best_log10p    = pmax(log10p_noBatch, log10p_Batch)
  )

pval_cor_pearson  <- cor(pval_compare$log10p_noBatch, pval_compare$log10p_Batch, method = "pearson")
pval_cor_spearman <- cor(pval_compare$log10p_noBatch, pval_compare$log10p_Batch, method = "spearman")

max_abs_diff <- max(pval_compare$abs_diff, na.rm = TRUE)
max_diff_snp <- pval_compare %>% filter(abs_diff == max_abs_diff)
diff_range   <- range(pval_compare$log10p_diff, na.rm = TRUE)

cat("==================== P-VALUE CORRELATION (per-SNP) ====================\n")
cat(sprintf("N SNPs compared:                %d\n", nrow(pval_compare)))
cat(sprintf("Pearson correlation (-log10 p): %.4f\n", pval_cor_pearson))
cat(sprintf("Spearman correlation (-log10 p):%.4f\n", pval_cor_spearman))
cat(sprintf("Max |change| in -log10(p):      %.4f\n", max_abs_diff))
cat(sprintf("Range of change (noBatch->Batch):[%.4f, %.4f]\n", diff_range[1], diff_range[2]))
cat("=========================================================================\n\n")

print(max_diff_snp)
write.csv(pval_compare, "pval_comparison_per_snp.csv", row.names = FALSE)

#####################################################################
# 4. Significance thresholds (suggestive = 1e-6, genome-wide = 4e-7)
#####################################################################

sig_threshold         <- 4e-7   # genome-wide
suggestive_threshold  <- 1e-6

pval_compare <- pval_compare %>%
  mutate(
    is_suggestive = best_log10p > -log10(suggestive_threshold),
    is_genomewide = best_log10p > -log10(sig_threshold)
  )

#####################################################################
# 5. Figure: magnitude of p-value shift vs. SNP significance
#    (are the biggest-changing SNPs also the most significant ones?)
#####################################################################

sig_vs_diff_plot <- ggplot(pval_compare, aes(x = best_log10p, y = abs_diff)) +
  geom_point(alpha = 0.15, size = 0.6, color = "grey40") +
  geom_point(data = filter(pval_compare, is_suggestive),
             aes(x = best_log10p, y = abs_diff), color = "steelblue", alpha = 0.6, size = 1) +
  geom_vline(xintercept = -log10(suggestive_threshold), linetype = "dashed", color = "orange") +
  geom_vline(xintercept = -log10(sig_threshold), linetype = "dashed", color = "red") +
  labs(
    x = expression("Most significant "*-log[10](P)*" across both models"),
    y = expression("Absolute change in "*-log[10](P)),
    title = "Magnitude of p-value shift vs. SNP significance",
    subtitle = sprintf("Orange line = suggestive (%.0e), red line = genome-wide (%.0e)",
                        suggestive_threshold, sig_threshold)
  ) +
  theme_minimal(base_size = 13)

print(sig_vs_diff_plot)
ggsave("sig_vs_diff_scatter.png", sig_vs_diff_plot, width = 7, height = 6, dpi = 300)

#####################################################################
# 6. Numeric summaries used in the supplementary text
#####################################################################

# Where does the single largest-shift SNP rank in significance?
rank_of_max_diff_snp <- pval_compare %>%
  arrange(desc(abs_diff)) %>%
  slice(1)

# SNPs reaching suggestive/genome-wide significance in either model
sig_snps_summary <- pval_compare %>%
  filter(is_suggestive) %>%
  summarise(
    n_snps        = n(),
    max_abs_diff  = max(abs_diff, na.rm = TRUE),
    mean_abs_diff = mean(abs_diff, na.rm = TRUE)
  )

# SNPs with at least a one-order-of-magnitude shift (abs_diff >= 1)
large_diff_summary_10x <- pval_compare %>%
  filter(abs_diff >= 1) %>%
  summarise(
    n_snps               = n(),
    min_abs_diff         = min(abs_diff, na.rm = TRUE),
    max_abs_diff         = max(abs_diff, na.rm = TRUE),
    mean_abs_diff        = mean(abs_diff, na.rm = TRUE),
    n_suggestive_or_above = sum(is_suggestive),
    n_genomewide         = sum(is_genomewide),
    min_p_value          = min(10^(-best_log10p), na.rm = TRUE)
  )

cat("==================== TOP SHIFTED SNP ====================\n")
print(rank_of_max_diff_snp)

cat("\n==================== SIGNIFICANT SNPs (suggestive/genome-wide) ====================\n")
print(sig_snps_summary)

cat("\n==================== SNPs WITH >=1 ORDER-OF-MAGNITUDE SHIFT ====================\n")
print(large_diff_summary_10x)