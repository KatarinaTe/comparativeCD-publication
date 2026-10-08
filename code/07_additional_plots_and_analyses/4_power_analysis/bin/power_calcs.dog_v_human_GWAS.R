#!/usr/bin/env Rscript
# power_calcs.dog_v_human_GWAS.R -- statistical power comparison, dogCD GWAS vs. human OCD GWAS
# (Strom et al. 2025 meta-analysis). DogCD achieves power comparable to human OCD despite a far
# smaller sample size, because canine breed structure inflates per-locus variance explained.
#
# Source: Elinor's power_calcs.dog_v_human_GWAS.R, provided as a self-contained bundle (script +
# 2 input TSVs + pre-computed outputs + a figure-legend/methods docx) under
# data/Elinor_simple_power_analysis/. Moved here 2026-09-24, following 2_gwas_catalog_overlap's
# precedent for a frozen-input, no-pipeline-prerequisite substage. Only change from the original:
# the two hardcoded input filenames below became CLI args -- no other logic touched. Original
# bundle (including the figure-legend/methods docx and pre-computed outputs, used to verify this
# conversion reproduces them exactly) preserved at
# archive/elinor_power_analysis_originals/ for provenance.

library(tidyverse)
library(cowplot)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
N_CASE   <- 53660        # Strom et al. 2025 OCD meta-analysis
N_CTRL   <- 2044417
P_HUMAN  <- 5e-8         # human genome-wide significance threshold
P_DOG    <- 4e-7         # dogCD discovery threshold (extended canine LD)
TARGET_POWER <- 0.80

NEFF_HUMAN <- 4 / (1 / N_CASE + 1 / N_CTRL)

# File paths -- passed in as CLI args (was hardcoded relative filenames in the original)
args <- commandArgs(trailingOnly = TRUE)
DOG_FILE   <- args[1]
HUMAN_FILE <- args[2]

# ---------------------------------------------------------------------------
# Load dogCD top hits
# columns include: Item number | SNP | P | A1/A2 | Freq (A2) | Beta (A2) | SE (A2) | Total (n)
# (plus clumping/region metadata columns, which are ignored here)
# ---------------------------------------------------------------------------
dog_raw <- read_tsv(DOG_FILE, show_col_types = FALSE)

dog <- dog_raw %>%
  rename(
    item  = `Item number`,
    snp   = SNP,
    p     = P,
    a1_a2 = `A1/A2`,
    freq  = `Freq (A2)`,
    beta  = `Beta (A2)`,
    se    = `SE (A2)`,
    n     = `Total (n)`
  ) %>%
  # keep every independently clumped region per trait (not just each trait's
  # single best locus) -- drop rows with no effect estimate (e.g. a region
  # with no lead SNP call)
  filter(!is.na(freq), !is.na(beta)) %>%
  distinct(item, snp, .keep_all = TRUE) %>%
  mutate(
    # R2 from the test statistic implied by the reported p-value (chi2 =
    # qchisq(1-p, 1)), not from 2pq*beta^2. This is model-agnostic -- it
    # doesn't assume beta is on a linear phenotype scale, so it's valid for
    # both the linear (MLMA-LOCO) and ordinal (POLMM) traits alike, and it
    # avoids relying on beta/SE, which are rounded in this file and don't
    # always exactly reproduce the reported p-value.
    chi2 = qchisq(1 - p, df = 1),
    r2 = chi2 / (chi2 + n),
    passes_threshold = p <= P_DOG
  )

# ---------------------------------------------------------------------------
# Load human GWAS table (Strom et al. 2025, Table 1: 30 genome-wide-significant loci)
# EDIT THIS MAPPING if your Strom.Table1.tsv uses different column headers.
# Expected fields: SNP, P value, OR, FRQCA (freq in cases), FRQCO (freq in controls)
# ---------------------------------------------------------------------------
human_raw <- read_tsv(HUMAN_FILE, show_col_types = FALSE)

required_human_cols <- c("SNP", "P value", "OR", "FRQCA", "FRQCO")
missing_cols <- setdiff(required_human_cols, names(human_raw))
if (length(missing_cols) > 0) {
  stop(
    "Strom.Table1.tsv is missing expected column(s): ", paste(missing_cols, collapse = ", "),
    "\nUpdate the `rename()` block below to match your file's actual headers: ",
    paste(names(human_raw), collapse = ", ")
  )
}

human <- human_raw %>%
  rename(
    snp   = SNP,
    p     = `P value`,
    or    = OR,
    frqca = FRQCA,
    frqco = FRQCO
  ) %>%
  mutate(
    # Same test-statistic-derived R2 as dogCD, using Neff in place of N.
    chi2 = qchisq(1 - p, df = 1),
    r2 = chi2 / (chi2 + NEFF_HUMAN)
  )

# ---------------------------------------------------------------------------
# Power functions (non-central chi-square, 1 df)
# ---------------------------------------------------------------------------
power_at <- function(n, p_thresh, r2) {
  crit <- qchisq(1 - p_thresh, df = 1)
  ncp  <- n * r2 / (1 - r2)
  1 - pchisq(crit, df = 1, ncp = ncp)
}

r2_for_power <- function(n, p_thresh, target = TARGET_POWER) {
  r2_grid <- 10^seq(-6, -0.3, length.out = 20000)
  pw <- power_at(n, p_thresh, r2_grid)
  r2_grid[which(pw >= target)[1]]
}

# ---------------------------------------------------------------------------
# Per-trait / per-locus summary stats (for the text)
# ---------------------------------------------------------------------------
dog_summary <- dog %>%
  rowwise() %>%
  mutate(
    power_at_hit   = power_at(n, P_DOG, r2),
    r2_for_80power = r2_for_power(n, P_DOG)
  ) %>%
  ungroup() %>%
  select(item, snp, freq, beta, p, n, r2, passes_threshold, power_at_hit, r2_for_80power) %>%
  arrange(p)

human_summary <- human %>%
  rowwise() %>%
  mutate(power_at_hit = power_at(NEFF_HUMAN, P_HUMAN, r2)) %>%
  ungroup() %>%
  select(snp, frqco, or, p, r2, power_at_hit) %>%
  arrange(p)

r2_80_human <- r2_for_power(NEFF_HUMAN, P_HUMAN)

cat("Neff (human OCD, case-control) =", format(NEFF_HUMAN, big.mark = ","), "\n")
cat("R2 needed for 80% power, human OCD design:", signif(r2_80_human * 100, 3), "%\n")
cat("R2 needed for 80% power, dogCD designs (range):",
    signif(min(dog_summary$r2_for_80power) * 100, 3), "-",
    signif(max(dog_summary$r2_for_80power) * 100, 3), "%\n")
cat("Ratio (median dogCD / human):",
    signif(median(dog_summary$r2_for_80power) / r2_80_human, 1), "x\n\n")

print(dog_summary, n = Inf)
print(human_summary, n = Inf)

write_csv(dog_summary,   "dogCD_power_summary.csv")
write_csv(human_summary, "human_OCD_power_summary.csv")

# ---------------------------------------------------------------------------
# Build power curves
# dogCD is collapsed to a single representative curve using the minimum N
# across the 17 traits (most conservative choice; N ranges 2,322-2,555).
# Per-trait points still use each trait's own actual N/R2, not this minimum.
# ---------------------------------------------------------------------------
r2_grid <- 10^seq(-5, -0.3, length.out = 400)

N_DOG_REPR <- min(dog$n)

dog_curve <- tibble(r2 = r2_grid) %>%
  mutate(power = power_at(N_DOG_REPR, P_DOG, r2))

human_curve <- tibble(r2 = r2_grid) %>%
  mutate(power = power_at(NEFF_HUMAN, P_HUMAN, r2))

# ---------------------------------------------------------------------------
# Plot: one collapsed dogCD line + one human line, with an in-plot legend.
# Points are only plotted for loci that actually cross their own study's
# genome-wide significance threshold; dog and human get distinct filled
# shapes (circle vs. diamond), mapped through a shared "species" variable so
# color/linewidth/shape combine into a single legend.
# ---------------------------------------------------------------------------
DOG_COLOR   <- "#B2182B"  # sophisticated red
HUMAN_COLOR <- "#4D4D4D"  # dark grey, distinct from red against white
DOG_LABEL   <- "DogCD"
HUMAN_LABEL <- "Human OCD"

dog_curve      <- dog_curve      %>% mutate(species = DOG_LABEL)
dog_hits_sig   <- dog_summary    %>% filter(passes_threshold) %>% mutate(species = DOG_LABEL)
human_curve    <- human_curve    %>% mutate(species = HUMAN_LABEL)
human_hits_sig <- human_summary  %>% filter(p <= P_HUMAN)     %>% mutate(species = HUMAN_LABEL)

species_colors    <- c(DogCD = DOG_COLOR, `Human OCD` = HUMAN_COLOR)
species_shapes    <- c(DogCD = 16,        `Human OCD` = 18)          # filled circle / filled diamond
species_linetypes <- c(DogCD = "solid",   `Human OCD` = "solid")     # both solid; widths differ below
species_linewidth <- c(DogCD = 1.1,       `Human OCD` = 0.5)         # human thinner than dog

p_plot <- ggplot() +
  geom_line(
    data = bind_rows(dog_curve, human_curve),
    aes(x = r2 * 100, y = power, color = species, linetype = species, linewidth = species)
  ) +
  geom_point(
    data = bind_rows(dog_hits_sig, human_hits_sig),
    aes(x = r2 * 100, y = power_at_hit, color = species, shape = species),
    size = 3.6, alpha = 0.25
  ) +
  scale_x_log10(labels = scales::label_number(accuracy = 0.001, drop0trailing = TRUE)) +
  scale_color_manual(values = species_colors, name = NULL) +
  scale_shape_manual(values = species_shapes, name = NULL) +
  scale_linetype_manual(values = species_linetypes, name = NULL) +
  scale_linewidth_manual(values = species_linewidth, name = NULL) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(
    x = "Variance explained by locus (%)",
    y = "Power to reach discovery threshold",
    #title = "Statistical power: dogCD GWAS vs. human OCD GWAS",
    #subtitle = sprintf(
    #  "DogCD (N=%.0f, P<%.0e) vs. human OCD (Neff=%s, P<%.0e)",
    #  N_DOG_REPR, P_DOG, format(round(NEFF_HUMAN), big.mark = ","), P_HUMAN
    #)
  ) +
  theme_cowplot(font_size = 8) +
  theme(
    legend.position = c(0.03, 0.97),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "grey50"),
    legend.box.background = element_rect(fill = "white", color = "grey50"),
    legend.margin = margin(4, 6, 4, 6)
  )

ggsave("dogCD_vs_humanOCD.GWAS_power_comparison.png", p_plot, width = 4.25, height = 3, bg="white",dpi = 300)
ggsave("dogCD_vs_humanOCD.GWAS_power_comparison.pdf", p_plot, width = 4.25, height = 3)
p_plot
