#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
supp_tables_file <- args[1]
pdf_out          <- args[2]

# Human OCD SNP-heritability reference lines (ref. 1), drawn as dashed line +/- SE band.
# Per KatarinaTe 2026-10-06: all ascertainment types 6.7% +/- 0.3%; clinical ascertainment
# 16.4% +/- 1.5%.
human_h2 <- tibble(
  label = c("human\n(clinical\nascertain-\nment)", "human\n(all ascertain-\nment types)"),
  h2    = c(0.164, 0.067),
  se    = c(0.015, 0.003)
)

# Per KatarinaTe (PR #4 review comment, 2026-09-14): cite the published Supplementary Tables
# directly rather than the private Google Sheet, which won't be available to others.
heritabilities <- read_excel(supp_tables_file, sheet = "ST4", skip = 1, na = "NA")
heritabilities <- heritabilities %>% filter(!is.na(`Item number`))
question_key <- read_excel(supp_tables_file, sheet = "ST1", skip = 1, na = "")
question_key <- question_key %>% filter(!is.na(`Survey question`))
question_key$`Item number` <- as.numeric(question_key$`Item number`)

heritabilities <- heritabilities %>% select(`Item number`, `GREML LDS constrained heritability`, `GREML LDS noConstraint heritability`)
names(heritabilities) <- c("phenotype", "LDS_constrained", "LDS_no_constraint")
heritabilities$pheno_numeric <- as.numeric(gsub("CCDF|item", "", heritabilities$phenotype))

# ST1's own column is "Highest loading factor" (numeric 1/2/3), not "Factor" -- fixed 2026-09-22
# after the join started silently producing no `Factor` column at all against the current
# NATURE_SupplementaryTables.xlsx (a real schema change: this column used to be named "Factor"
# in an earlier version of ST1, per this script's own original assumption).
items <- heritabilities %>%
  filter(grepl("^item", phenotype), !is.na(LDS_constrained)) %>%
  left_join(question_key, by = c("pheno_numeric" = "Item number")) %>%
  mutate(Factor = paste("Factor", `Highest loading factor`), heritability = LDS_constrained)

# Factor-level heritabilities (stars): GREML LDS constrained, except Factor 2, whose LD-constrained
# model failed to converge (NA in ST4) -- it uses GREML LDS noConstraint (per KatarinaTe 2026-10-06).
factors <- heritabilities %>%
  filter(grepl("^CCDF", phenotype)) %>%
  mutate(Factor = paste("Factor", pheno_numeric),
         heritability = if_else(phenotype == "CCDF2", LDS_no_constraint, LDS_constrained))
stopifnot(!anyNA(factors$heritability))

# The two most heritable factor-3 items are labeled
labeled_items <- items %>%
  filter(`Highest loading factor` == 3) %>%
  slice_max(heritability, n = 2) %>%
  mutate(label = paste0("Q", pheno_numeric, "\n",
                        str_to_sentence(str_remove(str_remove(`Short name`, "^time spent "), " or circling$"))))

box_fill  <- "#3d9dbf"
box_line  <- "#286677"
text_grey <- "grey21"

p <- ggplot(items, aes(x = Factor, y = heritability)) +
  # human reference bands/lines
  geom_rect(data = human_h2, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = h2 - se, ymax = h2 + se),
            fill = "grey60", alpha = 0.25) +
  geom_hline(data = human_h2, aes(yintercept = h2), linetype = "dashed", linewidth = 0.4) +
  geom_text(data = human_h2, aes(x = 4.05, y = h2 + 0.006, label = label), hjust = 0, vjust = 1,
            size = 5 / .pt, lineheight = 0.9, colour = text_grey) +
  # per-item boxplots + points
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = box_fill, alpha = 0.5,
               colour = box_line, linewidth = 0.4) +
  stat_summary(fun = median, geom = "crossbar", width = 0.55, colour = "#122e35", linewidth = 0.3) +
  geom_point(size = 1.2, alpha = 0.7, colour = "grey25") +
  # factor-level heritability as star
  geom_point(data = factors, colour = "#36234a", shape = 8, size = 2.5, stroke = 0.4) +
  # labels for the two most heritable factor-3 items
  geom_segment(data = labeled_items, aes(x = 3.08, xend = 3.2, yend = heritability),
               linewidth = 0.25, colour = text_grey) +
  geom_text(data = labeled_items, aes(x = 3.25, label = label), hjust = 0, vjust = 0.75,
            size = 5 / .pt, lineheight = 0.9, colour = text_grey) +
  scale_x_discrete() +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3)) +
  coord_cartesian(xlim = c(0.6, 3.4), clip = "off") +
  labs(x = NULL, y = expression("heritability (" * h[SNP]^2 * ")")) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 5),
        axis.line = element_line(colour = "black", linewidth = 0.3),
        axis.ticks = element_line(colour = "black", linewidth = 0.3),
        plot.margin = margin(4, 60, 2, 2))

ggsave(pdf_out, p, width = 2.25, height = 2.4)
