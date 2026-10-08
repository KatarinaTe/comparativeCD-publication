#!/usr/bin/env Rscript
# plot_netcoloc_zscores.R -- Fig. 1f: dogCD (NPS_d) vs. human OCD (NPS_h) network-propagation
# z-score scatter plot, colored by seed-gene membership and conserved-network status, plus a
# diagnostic boxplot of combined z-scores by TOGA orthology class.
#
# Source: 05_network/1_NetColoc/plot_netcoloc_zscores.R (logic unmodified from the original,
# only file I/O parameterized -- see below). Labeled "Fig 1e" in the original script/notebook
# (1.1_OCD_dogCD_NetColoc_analysis_260521.ipynb cell 25: "make scatter plot for Figure 1e"),
# before the manuscript's panel lettering shifted to the current "Fig 1f".
#
# Two changes from the original:
#   1. Reads 1_NetColoc's own COMBINE_ZSCORES output (*_CCD_zcomb_z12.tsv: gene, NPS_r [dogCD],
#      NPS_h [human trait], NPS_hr [combined]) directly, instead of recomputing the same
#      full_join(ocd, ccd) + D1_z*D2_z multiplication the original did from raw per-trait
#      z-score files -- COMBINE_ZSCORES already produces byte-identical values (confirmed
#      against the deposited data/05_network/ocd_ccd_zcomb_z12_251022.txt), so this avoids
#      duplicating that computation. NPS_h/NPS_r renamed back to D1_z/D2_z below purely so the
#      rest of the original logic (column names, thresholds) needs no further changes.
#   2. `library(eulerr)` dropped -- loaded in the original but never actually called anywhere in
#      this script (likely a leftover from copying Fig2_network-overlap.R, which does use it).
#      `size=` on geom_vline/geom_hline/geom_segment changed to `linewidth=` (ggplot2 >=3.4
#      deprecation; same visual result).
#
# setwd() and every hardcoded absolute path from the original removed -- all inputs/outputs are
# now CLI args (paths staged by Nextflow; see process_definitions.nf).

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
zcomb_file         <- args[1]  # <trait>_CCD_zcomb_z12.tsv: gene, NPS_r (dogCD), NPS_h (OCD), NPS_hr
toga_file          <- args[2]
human_seeds_file   <- args[3]
dog_seeds_file     <- args[4]
combined_table_out <- args[5]
boxplot_out        <- args[6]
source_data_out    <- args[7]
scatter_out        <- args[8]

toga <- read.csv(toga_file, sep = "\t", header = TRUE)

zcomb <- read.csv(zcomb_file, sep = "\t", header = TRUE)
ocd_ccd <- data.frame(gene = zcomb$gene, D1_z = zcomb$NPS_h, D2_z = zcomb$NPS_r, zcomb = zcomb$NPS_hr)

write.table(ocd_ccd, file = combined_table_out, sep = "\t", row.names = FALSE, quote = FALSE)

#### check orthologs ####
genestoga <- full_join(ocd_ccd, toga, by = "gene")
genestoga2 <- na.omit(genestoga)

# Count the number of genes in each orthology class
class_counts <- genestoga2 %>%
  group_by(orthology_class) %>%
  summarize(n = n()) %>%
  mutate(label = paste0(orthology_class, "\n(n = ", n, ")"))

# Add the CONSERVED_NETWORK column to genestoga2
genestoga2$CONSERVED_NETWORK <- ifelse(genestoga2$D2_z > 1.5 & genestoga2$D1_z > 1.5 & genestoga2$zcomb > 3, 1, 0)

# Diagnostic plot: distribution of combined z-scores by orthology class (feeds nothing
# downstream -- kept as a reproduced leaf output, same as other diagnostic plots in this repo)
p_box <- ggplot(genestoga2, aes(x = orthology_class, y = zcomb)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = factor(CONSERVED_NETWORK)),
              width = 0.2, alpha = 0.5) +
  geom_text_repel(data = subset(genestoga2, CONSERVED_NETWORK == 1),
                  aes(label = gene),
                  size = 3,
                  box.padding = 0.5,
                  point.padding = 0.2,
                  force = 2,
                  segment.color = NA) +
  scale_x_discrete(labels = class_counts$label) +
  scale_color_manual(values = c("0" = "grey", "1" = "red"),
                     name = "Conserved Network",
                     labels = c("No", "Yes")) +
  theme_bw() +
  labs(x = "Orthology Class",
       y = "Combined Z-score",
       title = "Distribution of Combined Z-scores by Orthology Class") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
ggsave(boxplot_out, p_box)

#### Fig 1f scatter plot ####

# Hardcoded gene lists to label on the plot -- unchanged from the original (a fixed, specific
# subset of the conserved-network genes chosen for readability, not a general input).
dog_seeds_network <- c("NTRK2", "PROSER3", "GABRG1", "CDH8", "GRM7", "DACH1", "KLHL1", "KATNIP")
human_seeds_network <- c("ARIH2", "P4HTM", "KLHDC8B", "LSAMP", "CCHCR1", "LRFN5")

dog_seeds <- scan(dog_seeds_file, what = "", sep = "\n")
human_seeds <- scan(human_seeds_file, what = "", sep = "\n")

# Define coloring logic with combined condition for orange
genestoga2$color <- "grey"
genestoga2$color[genestoga2$gene %in% human_seeds] <- "cornflowerblue"
genestoga2$color[genestoga2$gene %in% dog_seeds] <- "coral"
genestoga2$color[
  genestoga2$zcomb > 3 & genestoga2$D1_z > 1.5 & genestoga2$D2_z > 1.5
] <- "orange"  # orange if all these conditions met (overrides others)

# Create annotation column for seed genes (unchanged)
genestoga2$label <- ""
genestoga2$label[genestoga2$gene %in% human_seeds_network] <- genestoga2$gene[genestoga2$gene %in% human_seeds_network]
genestoga2$label[genestoga2$gene %in% dog_seeds_network] <- genestoga2$gene[genestoga2$gene %in% dog_seeds_network]

write.table(genestoga2, file = source_data_out, sep = "\t", row.names = FALSE, quote = FALSE)

# Plot
p_scatter <- ggplot(genestoga2, aes(x = D1_z, y = D2_z)) +
  geom_point(aes(color = color), size = 3, alpha = 0.8) +
  scale_color_identity() +
  geom_text(aes(label = label), hjust = -0.1, vjust = 0.3, size = 3) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 1.5, xend = 1.5, y = 1.5, yend = Inf),
               linetype = "dotted", color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 1.5, xend = Inf, y = 1.5, yend = 1.5),
               linetype = "dotted", color = "black", linewidth = 0.7) +
  theme_minimal()
ggsave(scatter_out, p_scatter)
