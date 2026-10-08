#!/usr/bin/env Rscript
# plot_network_overlap.R -- Fig. 2b: unique/shared genes between the OCD/dogCD cross-species
# network and the human-only Depression/Schizophrenia networks (stacked bar + Euler diagram).
#
# Source: 05_network/4_Plotting/Fig2_network-overlap.R (logic unmodified, only file I/O
# parameterized -- see below).
#
# Panel A moved out 2026-09-27 into its own process (BUILD_NETWORK_GENE_COUNTS + PLOT_FIG2A,
# see bin/remake_fig2A.R) -- its 7 (captured, total) pairs were hardcoded literals here, and one
# was found genuinely wrong: Depression/dogCD's "total" was coded as 387 assuming the same
# z_comb=3 threshold as OCD/dogCD, but its real source notebook
# (DEP_CCD_Oct25_NetColoc_analysis_251022.ipynb) explicitly used z_comb=4 ("default = 3") --
# 387 is correct, but only at that different threshold, and the hand-verification process that
# caught the earlier "Depression human" copy-paste bug (2026-09-18) didn't catch this one. Fig.
# 2a is now fully live-computed from the same deposited z-score/seed-gene/systems-map files this
# script already uses for Panel B, rather than another hand-verified snapshot -- see
# REPRODUCIBILITY_AUDIT.md's 2026-09-27 entries for the full investigation.
#
# Panel B (stacked bar + Euler diagram) IS live-computed from the real systems-map gene lists.
# Fixed 2026-09-22: previously read dep_ccd_/sch_ccd_systemsmap_genes.txt (the CROSS-SPECIES
# Depression/Schizophrenia files) here by mistake -- the manuscript caption is explicit that this
# panel compares the cross-species OCD/dogCD network against the HUMAN-ONLY Depression/
# Schizophrenia gene sets, not their own cross-species counterparts. That bug silently produced
# the wrong overlap numbers (74/156/233 unique instead of 106/420/511, etc.) despite
# data/Manuscript_draft/REVIEW_AND_SUGGESTIONS.md's 2026-09-18 note claiming this panel had
# already been verified end-to-end. Re-confirmed 2026-09-22 with the corrected file references:
# all 7 Venn numbers (420/106/511 unique, 28/62/14 pairwise, 17 triple) now match the current
# published Fig. 2b image exactly.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(eulerr)
})

args <- commandArgs(trailingOnly = TRUE)
ocd_ccd_genes_file <- args[1]  # ocd_ccd_systemsmap_genes.txt -- cross-species OCD/dogCD
dep_genes_file      <- args[2]  # dep_systemsmap_genes.txt -- human-only Depression
scz_genes_file       <- args[3]  # sch_systemsmap_genes.txt -- human-only Schizophrenia
plot_out            <- args[4]

#### Panel B: live from real systems-map gene lists ####
ocd_genes <- readLines(ocd_ccd_genes_file)
dep_genes <- readLines(dep_genes_file)
scz_genes <- readLines(scz_genes_file)

ocd_tot <- length(ocd_genes); dep_tot <- length(dep_genes); scz_tot <- length(scz_genes)

shared_all_set <- Reduce(intersect, list(ocd_genes, dep_genes, scz_genes))
ocd_dep_set <- setdiff(intersect(ocd_genes, dep_genes), shared_all_set)
ocd_scz_set <- setdiff(intersect(ocd_genes, scz_genes), shared_all_set)
dep_scz_set <- setdiff(intersect(dep_genes, scz_genes), shared_all_set)

shared_all <- length(shared_all_set)
ocd_dep <- length(ocd_dep_set); ocd_scz <- length(ocd_scz_set); dep_scz <- length(dep_scz_set)

ocd_u <- length(setdiff(ocd_genes, union(dep_genes, scz_genes)))
dep_u <- length(setdiff(dep_genes, union(ocd_genes, scz_genes)))
scz_u <- length(setdiff(scz_genes, union(ocd_genes, dep_genes)))

tots3 <- c(ocd_tot, dep_tot, scz_tot)

df_b <- tibble(
  network  = rep(c("OCD/dogCD","Depression/dogCD","SCZ/dogCD"), 3),
  category = rep(c("Unique to network","Pairwise shared","Shared all three"), each=3),
  pct = c(
    round(c(ocd_u,dep_u,scz_u)/tots3*100, 1),
    round(c(ocd_dep+ocd_scz, ocd_dep+dep_scz, ocd_scz+dep_scz)/tots3*100, 1),
    round(shared_all/tots3*100, 1)
  )
) %>%
  mutate(
    category = factor(category, levels = c("Unique to network","Pairwise shared","Shared all three")),
    network  = factor(network,  levels = c("OCD/dogCD","Depression/dogCD","SCZ/dogCD")),
    n        = round(pct/100 * rep(tots3, 3))
  )

pB_bar <- ggplot(df_b, aes(x = network, y = pct, fill = category)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(round(pct), "%\n(", n, ")"),
                color = category),
            position = position_stack(vjust = 0.5), size = 2.8, lineheight = 1.1) +
  scale_fill_manual(values = c(
    "Unique to network" = "#534AB7",
    "Pairwise shared"   = "#AFA9EC",
    "Shared all three"  = "#EEEDFE")) +
  scale_color_manual(values = c(
    "Unique to network" = "#ffffff",
    "Pairwise shared"   = "#26215C",
    "Shared all three"  = "#7F77DD")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = NULL, y = "% of network genes", fill = NULL,
       title = "b   Unique and shared genes") +
  guides(color = "none") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        legend.key.size = unit(0.4, "cm"),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 11, face = "bold"))

#### Panel B: Euler diagram ####
fit <- euler(c(
  "OCD"         = ocd_u,
  "Dep"         = dep_u,
  "SCZ"         = scz_u,
  "OCD&Dep"     = ocd_dep,
  "OCD&SCZ"     = ocd_scz,
  "Dep&SCZ"     = dep_scz,
  "OCD&Dep&SCZ" = shared_all
))

pB_venn <- plot(fit,
  quantities = list(fontsize = 8, col = c("#534AB7","#0F6E56","#993C1D","grey30","grey30","grey30","#26215C")),
  labels     = list(labels = c("OCD/dogCD","Depression/dogCD","SCZ/dogCD"), fontsize = 8),
  fills      = list(fill = c("#EEEDFE","#E1F5EE","#FAECE7"), alpha = 0.65),
  edges      = list(col  = c("#534AB7","#0F6E56","#993C1D"), lwd = 1.2)
)

#### Combine ####
pB_bar / wrap_elements(pB_venn)

ggsave(plot_out, width = 5.5, height = 5.5, dpi = 300)
