library(tidyverse)

# Counts of brain and other tissue CREs from EPIC data
N_brain <- 126410 
N_other <- 238469 

# Analysis of overlap between clump +- 100Kb regions and EPIC dog CREs - do we see significant enrichment of active brain elements
# in these regions?


# 1. Set up the contingency table for testing significance
# O_brain  = Overlap with GWAS (Brain)
O_brain <- 408   
    
# O_other = Overlap with GWAS (Other)
O_other <- 655       
  
# Calculate the "No Overlap" groups
no_O_brain <- N_brain - O_brain
no_O_other <- N_other - O_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
contingency_table <- matrix(
  c(O_brain, no_O_brain,
    O_other, no_O_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
fisher_results <- fisher.test(contingency_table, alternative = "greater")

# 4. View the results
print(fisher_results)

# Answer: Yes! but not a huge effect: odds ratio = 1.17564, p = 0.0059

## Compare this to the SIZE regions

size_O_brain <- 325   
size_O_other <- 736       

# Calculate the "No Overlap" groups
size_no_O_brain <- N_brain - size_O_brain
size_no_O_other <- N_other - size_O_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
size_contingency_table <- matrix(
  c(size_O_brain, size_no_O_brain,
    size_O_other, size_no_O_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
size_fisher_results <- fisher.test(size_contingency_table, alternative = "greater")

# 4. View the results
print(size_fisher_results)

# Answer: Not significant! p = 0.9975, OR = 0.83

## Ok, now test for the GWS plus suggestive SNPs clumps

# First, for dogCD
# 1. Set up the contingency table for testing significance
dogCD_GWS_sugg_O_brain <- 698  
dogCD_GWS_sugg_O_other <- 1300

# Calculate the "No Overlap" groups
dogCD_GWS_sugg_no_O_brain <- N_brain - dogCD_GWS_sugg_O_brain
dogCD_GWS_sugg_no_O_other <- N_other - dogCD_GWS_sugg_O_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
dogCD_GWS_sugg_contingency_table <- matrix(
  c(dogCD_GWS_sugg_O_brain, dogCD_GWS_sugg_no_O_brain,
    dogCD_GWS_sugg_O_other, dogCD_GWS_sugg_no_O_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
dogCD_GWS_sugg_fisher_results <- fisher.test(dogCD_GWS_sugg_contingency_table, alternative = "greater")

# 4. View the results
print(dogCD_GWS_sugg_fisher_results)

# Answer:

## Compare this to the SIZE regions

size_GWS_sugg_O_brain <- 430   
size_GWS_sugg_O_other <- 996       


# Calculate the "No Overlap" groups
size_GWS_sugg_no_O_brain <- N_brain - size_GWS_sugg_O_brain
size_GWS_sugg_no_O_other <- N_other - size_GWS_sugg_O_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
size_GWS_sugg_contingency_table <- matrix(
  c(size_GWS_sugg_O_brain, size_GWS_sugg_no_O_brain,
    size_GWS_sugg_O_other, size_GWS_sugg_no_O_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
size_GWS_sugg_fisher_results <- fisher.test(size_GWS_sugg_contingency_table, alternative = "greater")

# 4. View the results
print(size_GWS_sugg_fisher_results)

# Answer: Not significant! p = 0.9999, OR = 0.81

######################


## Ok, let's redo this but rather than counting elements (which are highly variable in size), we count the base pair overlap, i.e. how many bps of our
# GWAS regions is covered by CREs
# Need to count the total bps of our brain and other regions CREs:

brain_cre_bp <- 247096600
other_cre_bp <- 393644600

# Now let's look at the DogCD GWS clumps +/- 100KB

dogCD_GWS_O_bp_brain <- 799607  
dogCD_GWS_O_bp_other <- 1109740      

# Calculate the "No Overlap" groups
dogCD_GWS_no_O_bp_brain <- brain_cre_bp - dogCD_GWS_O_bp_brain
dogCD_GWS_no_O_bp_other <- other_cre_bp - dogCD_GWS_O_bp_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
dogCD_GWS_bp_contingency_table <- matrix(
  c(dogCD_GWS_O_bp_brain, dogCD_GWS_no_O_bp_brain,
    dogCD_GWS_O_bp_other, dogCD_GWS_no_O_bp_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
dogCD_GWS_bp_fisher_results <- fisher.test(dogCD_GWS_bp_contingency_table, alternative = "greater")

# 4. View the results
print(dogCD_GWS_bp_fisher_results)

# Answer: p < 2.2x10-16, OR = 1.15

# Now let's look at the SIZE GWS clumps +/- 100KB

SIZE_GWS_O_bp_brain <- 710315
SIZE_GWS_O_bp_other <- 1330316

# Calculate the "No Overlap" groups
SIZE_GWS_no_O_bp_brain <- brain_cre_bp - SIZE_GWS_O_bp_brain
SIZE_GWS_no_O_bp_other <- other_cre_bp - SIZE_GWS_O_bp_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
SIZE_GWS_bp_contingency_table <- matrix(
  c(SIZE_GWS_O_bp_brain, SIZE_GWS_no_O_bp_brain,
    SIZE_GWS_O_bp_other, SIZE_GWS_no_O_bp_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
SIZE_GWS_bp_fisher_results <- fisher.test(SIZE_GWS_bp_contingency_table, alternative = "greater")

# 4. View the results
print(SIZE_GWS_bp_fisher_results)

# Answer: p = 1, OR = 0.85

#######################
#Now adding in UU regions:::
## Ok, let's redo this but rather than counting elements (which are highly variable in size), we count the base pair overlap, i.e. how many bps of our
# GWAS regions is covered by CREs
# Need to count the total bps of our brain and other regions CREs:

brain_cre_bp <- 247096600

# Counts when including UU regions
brain_cre_bp <- 460303000

other_cre_bp <- 393644600

# Now let's look at the DogCD GWS clumps +/- 100KB

dogCD_GWS_O_bp_brain <- 799607

# If including UU regions
dogCD_GWS_O_bp_brain <- 1472848

dogCD_GWS_O_bp_other <- 1109740      

# Calculate the "No Overlap" groups
dogCD_GWS_no_O_bp_brain <- brain_cre_bp - dogCD_GWS_O_bp_brain
dogCD_GWS_no_O_bp_other <- other_cre_bp - dogCD_GWS_O_bp_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
dogCD_GWS_bp_contingency_table <- matrix(
  c(dogCD_GWS_O_bp_brain, dogCD_GWS_no_O_bp_brain,
    dogCD_GWS_O_bp_other, dogCD_GWS_no_O_bp_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
dogCD_GWS_bp_fisher_results <- fisher.test(dogCD_GWS_bp_contingency_table, alternative = "greater")

# 4. View the results
print(dogCD_GWS_bp_fisher_results)

# Answer: p < 2.2x10-16, OR = 1.15

# Now let's look at the SIZE GWS clumps +/- 100KB

SIZE_GWS_O_bp_brain <- 710315

#With uu regions
SIZE_GWS_O_bp_brain <- 1272105

SIZE_GWS_O_bp_other <- 1330316

# Calculate the "No Overlap" groups
SIZE_GWS_no_O_bp_brain <- brain_cre_bp - SIZE_GWS_O_bp_brain
SIZE_GWS_no_O_bp_other <- other_cre_bp - SIZE_GWS_O_bp_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
SIZE_GWS_bp_contingency_table <- matrix(
  c(SIZE_GWS_O_bp_brain, SIZE_GWS_no_O_bp_brain,
    SIZE_GWS_O_bp_other, SIZE_GWS_no_O_bp_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
SIZE_GWS_bp_fisher_results <- fisher.test(SIZE_GWS_bp_contingency_table, alternative = "greater")

# 4. View the results
print(SIZE_GWS_bp_fisher_results)

# Answer: p = 1, OR = 0.85





## How about if we reduce the regions to just the 'clumps' regions without the +-100Kb

# 1. Set up the contingency table for testing significance
# O_brain  = Overlap with GWAS (Brain)
# N_brain  = Total Brain elements
O_brain <- 98   
N_brain <- 126410     

# O_other = Overlap with GWAS (Other)
# N_other = Total Other elements
O_other <- 173       
N_other <- 238469    

# Calculate the "No Overlap" groups
no_O_brain <- N_brain - O_brain
no_O_other <- N_other - O_other

# 2. Construct the 2x2 Contingency Table
# byrow = TRUE ensures the data fills row-by-row mapping to our dimnames
contingency_table <- matrix(
  c(O_brain, no_O_brain,
    O_other, no_O_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Tissue = c("Brain", "Other"),
    GWAS_Overlap = c("Yes", "No")
  )
)

# 3. Run the Fisher's Exact Test
# We use alternative="greater" because you want to know if Brain > Other. 
# If you just wanted to know if they were *different*, you'd use "two.sided" (the default).
fisher_results <- fisher.test(contingency_table, alternative = "greater")

# 4. View the results
print(fisher_results)
# No longer significant, odds ration = 1.07, p = 0.32




## Tissue-specificity figure: does GWAS region overlap with cCREs depend on
## tissue (brain vs. all other tissues), separately for dogCD and body size?
## Panel A: raw overlap rate, Brain vs. Other, per trait (dumbbell)
## Panel B: odds ratio (Brain vs. Other), 95% CI, per trait (Fisher's exact test)

library(ggplot2)
library(dplyr)
library(patchwork)

## ---- Inputs (bp-level overlap of GWAS regions with cCREs) ----
## epicdog_tissue_bp_overlaps.txt is real bedtools output (code/06_cCRE/check_epicdog_bp_overlap.sh,
## run against the fetched EpicDog per-tissue BEDs and both GWAS region sets), not hardcoded literals.
## Fixed 2026-09-17: the values these used to be hardcoded to had silently dropped mammary-gland
## tissue from the "other" pool (a mammary_gland/mammary filename mismatch in the original
## check_CRE_overlap.sh's `cat` step) — see data/Manuscript_draft/REVIEW_AND_SUGGESTIONS.md.
epicdog_overlaps <- read.delim("~/Documents/GitHub/comparativeCD/data/06_cCRE/epicdog_tissue_bp_overlaps.txt")

brain_cre_bp <- epicdog_overlaps$Total_bp[epicdog_overlaps$Tissue_category == "brain"][1]
other_cre_bp <- epicdog_overlaps$Total_bp[epicdog_overlaps$Tissue_category == "other"][1]

dogCD_O_brain <- epicdog_overlaps$Overlap_bp[epicdog_overlaps$GWAS_set == "dogCD" & epicdog_overlaps$Tissue_category == "brain"]
dogCD_O_other <- epicdog_overlaps$Overlap_bp[epicdog_overlaps$GWAS_set == "dogCD" & epicdog_overlaps$Tissue_category == "other"]
SIZE_O_brain  <- epicdog_overlaps$Overlap_bp[epicdog_overlaps$GWAS_set == "SIZE"  & epicdog_overlaps$Tissue_category == "brain"]
SIZE_O_other  <- epicdog_overlaps$Overlap_bp[epicdog_overlaps$GWAS_set == "SIZE"  & epicdog_overlaps$Tissue_category == "other"]

dogCD_tab <- matrix(c(dogCD_O_brain, brain_cre_bp - dogCD_O_brain,
                      dogCD_O_other, other_cre_bp - dogCD_O_other),
                    nrow = 2, byrow = TRUE,
                    dimnames = list(Tissue = c("Brain","Other"), GWAS_Overlap = c("Yes","No")))

SIZE_tab <- matrix(c(SIZE_O_brain, brain_cre_bp - SIZE_O_brain,
                     SIZE_O_other, other_cre_bp - SIZE_O_other),
                   nrow = 2, byrow = TRUE,
                   dimnames = list(Tissue = c("Brain","Other"), GWAS_Overlap = c("Yes","No")))

dogCD_fit <- fisher.test(dogCD_tab)  # two-sided, for a plottable finite 95% CI
SIZE_fit  <- fisher.test(SIZE_tab)

## ---- Panel A: raw overlap rates (dumbbell, not a 0-baseline bar, since the
## rates are tightly clustered between ~0.28% and ~0.34% and a 0-anchored bar
## would make all four values look visually identical) ----
rates_df <- data.frame(
  trait  = factor(c("Body size","Body size","dogCD","dogCD"), levels = c("Body size","dogCD")),
  tissue = factor(c("Brain","Other","Brain","Other"), levels = c("Other","Brain")),
  rate   = c(SIZE_O_brain / brain_cre_bp, SIZE_O_other / other_cre_bp,
             dogCD_O_brain / brain_cre_bp, dogCD_O_other / other_cre_bp)
)

p_rates <- ggplot(rates_df, aes(x = rate, y = trait)) +
  geom_line(aes(group = trait), color = "grey55", linewidth = 0.6) +
  geom_point(aes(color = tissue), size = 3.6) +
  scale_color_manual(values = c("Brain" = "#2C5F7C", "Other" = "grey60"), name = NULL) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 0.01)) +
  labs(x = "Overlap rate (GWS region bp in cCRE bp)", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.margin = margin(b = -4, l = -2),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

## ---- Panel B: odds ratio, Brain vs. Other, per trait ----
or_df <- data.frame(
  trait   = factor(c("Body size","dogCD"), levels = c("Body size","dogCD")),
  OR      = c(unname(SIZE_fit$estimate), unname(dogCD_fit$estimate)),
  CI_low  = c(SIZE_fit$conf.int[1], dogCD_fit$conf.int[1]),
  CI_high = c(SIZE_fit$conf.int[2], dogCD_fit$conf.int[2])
) %>%
  mutate(label = sprintf("%.3f (%.3f\u2013%.3f)", OR, CI_low, CI_high))

p_or <- ggplot(or_df, aes(y = trait)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high, color = trait),
                orientation = "y", width = 0.35, linewidth = 0.9) +
  geom_point(aes(x = OR, color = trait), size = 2.0) +
  geom_text(aes(x = 1.18, label = label), hjust = 0, size = 2.95, family = "sans", color = "grey20") +
  scale_color_manual(values = c("dogCD" = "#D55E00", "Body size" = "grey40"), guide = "none") +
  scale_x_continuous(limits = c(0.78, 1.55), breaks = seq(0.8, 1.4, 0.2)) +
  labs(x = "Odds ratio (Brain vs. other tissues)", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    plot.margin = margin(t = 23),  # align with panel A under its legend
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

combined <- p_rates + p_or +
  plot_layout(widths = c(1, 1.25)) +
  plot_annotation(
    title = "dogCD, but not body-size, GWAS regions are enriched in brain cCREs",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 12, family = "sans"))
  )

ggsave("tissue_specificity_plot.pdf", combined, width = 9.0, height = 2.6, device = cairo_pdf)
ggsave("tissue_specificity_plot.png", combined, width = 9.0, height = 2.6, dpi = 300)

cat("dogCD OR/CI: ", dogCD_fit$estimate, dogCD_fit$conf.int, "\n")
cat("SIZE  OR/CI: ", SIZE_fit$estimate, SIZE_fit$conf.int, "\n")

## Tissue-specificity figure: does GWAS region overlap with cCREs depend on
## tissue (brain vs. all other tissues), separately for dogCD and body size?
## Panel A: raw overlap rate, Brain vs. Other, per trait (dumbbell)
## Panel B: odds ratio (Brain vs. Other), 95% CI, per trait (Fisher's exact test)

library(ggplot2)
library(dplyr)
library(patchwork)

## ---- Inputs (bp-level overlap of GWAS regions with cCREs) ----
## epicdog_tissue_bp_overlaps.txt is real bedtools output (code/06_cCRE/check_epicdog_bp_overlap.sh,
## run against the fetched EpicDog per-tissue BEDs and both GWAS region sets), not hardcoded literals.
## Fixed 2026-09-17: the values these used to be hardcoded to had silently dropped mammary-gland
## tissue from the "other" pool (a mammary_gland/mammary filename mismatch in the original
## check_CRE_overlap.sh's `cat` step) — see data/Manuscript_draft/REVIEW_AND_SUGGESTIONS.md.
epicdog_overlaps <- read.delim("~/Documents/GitHub/comparativeCD/data/06_cCRE/epicdog_tissue_bp_overlaps.txt")

brain_cre_bp <- epicdog_overlaps$Total_bp[epicdog_overlaps$Tissue_category == "brain"][1]
other_cre_bp <- epicdog_overlaps$Total_bp[epicdog_overlaps$Tissue_category == "other"][1]

dogCD_O_brain <- epicdog_overlaps$Overlap_bp[epicdog_overlaps$GWAS_set == "dogCD" & epicdog_overlaps$Tissue_category == "brain"]
dogCD_O_other <- epicdog_overlaps$Overlap_bp[epicdog_overlaps$GWAS_set == "dogCD" & epicdog_overlaps$Tissue_category == "other"]
SIZE_O_brain  <- epicdog_overlaps$Overlap_bp[epicdog_overlaps$GWAS_set == "SIZE"  & epicdog_overlaps$Tissue_category == "brain"]
SIZE_O_other  <- epicdog_overlaps$Overlap_bp[epicdog_overlaps$GWAS_set == "SIZE"  & epicdog_overlaps$Tissue_category == "other"]

dogCD_tab <- matrix(c(dogCD_O_brain, brain_cre_bp - dogCD_O_brain,
                      dogCD_O_other, other_cre_bp - dogCD_O_other),
                    nrow = 2, byrow = TRUE,
                    dimnames = list(Tissue = c("Brain","Other"), GWAS_Overlap = c("Yes","No")))

SIZE_tab <- matrix(c(SIZE_O_brain, brain_cre_bp - SIZE_O_brain,
                     SIZE_O_other, other_cre_bp - SIZE_O_other),
                   nrow = 2, byrow = TRUE,
                   dimnames = list(Tissue = c("Brain","Other"), GWAS_Overlap = c("Yes","No")))

dogCD_fit <- fisher.test(dogCD_tab)  # two-sided, for a plottable finite 95% CI
SIZE_fit  <- fisher.test(SIZE_tab)

## ---- Panel A: raw overlap rates (dumbbell, not a 0-baseline bar, since the
## rates are tightly clustered between ~0.28% and ~0.34% and a 0-anchored bar
## would make all four values look visually identical) ----
rates_df <- data.frame(
  trait  = factor(c("Body size","Body size","dogCD","dogCD"), levels = c("Body size","dogCD")),
  tissue = factor(c("Brain","Other tissues","Brain","Other tissues"), levels = c("Other tissues","Brain")),
  rate   = c(SIZE_O_brain / brain_cre_bp, SIZE_O_other / other_cre_bp,
             dogCD_O_brain / brain_cre_bp, dogCD_O_other / other_cre_bp)
)

p_rates <- ggplot(rates_df, aes(x = rate, y = trait, color = tissue)) +
  geom_point(position = position_dodge(width = 0.5), size = 3.8) +
  scale_color_manual(values = c("Brain" = "#2C5F7C", "Other tissues" = "grey55"), name = NULL) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 0.01)) +
  labs(x = "Overlap rate (GWS region bp in cCRE bp)", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.margin = margin(b = -4, l = -2),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

## ---- Panel B: odds ratio, Brain vs. Other, per trait ----
or_df <- data.frame(
  trait   = factor(c("Body size","dogCD"), levels = c("Body size","dogCD")),
  OR      = c(unname(SIZE_fit$estimate), unname(dogCD_fit$estimate)),
  CI_low  = c(SIZE_fit$conf.int[1], dogCD_fit$conf.int[1]),
  CI_high = c(SIZE_fit$conf.int[2], dogCD_fit$conf.int[2])
) %>%
  mutate(label = sprintf("%.3f (%.3f\u2013%.3f)", OR, CI_low, CI_high))

p_or <- ggplot(or_df, aes(y = trait)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high, color = trait),
                orientation = "y", width = 0.35, linewidth = 0.9) +
  geom_point(aes(x = OR, color = trait), size = 2.0) +
  geom_text(aes(x = 1.18, label = label), hjust = 0, size = 2.95, family = "sans", color = "grey20") +
  scale_color_manual(values = c("dogCD" = "#D55E00", "Body size" = "grey40"), guide = "none") +
  scale_x_continuous(limits = c(0.78, 1.55), breaks = seq(0.8, 1.4, 0.2)) +
  labs(x = "Odds ratio (Brain vs. other tissues)", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    plot.margin = margin(t = 23),  # align with panel A under its legend
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

combined <- p_rates + p_or +
  plot_layout(widths = c(1, 1.25)) +
  plot_annotation(
    title = "dogCD, but not body-size, GWAS regions are enriched in brain cCREs",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 12, family = "sans"))
  )

ggsave("tissue_specificity_plot.pdf", combined, width = 9.0, height = 2.6, device = cairo_pdf)
ggsave("tissue_specificity_plot.png", combined, width = 9.0, height = 2.6, dpi = 300)

cat("dogCD OR/CI: ", dogCD_fit$estimate, dogCD_fit$conf.int, "\n")
cat("SIZE  OR/CI: ", SIZE_fit$estimate, SIZE_fit$conf.int, "\n")

#### UU dataset ####
## Ok, now we know there is significant enrichment (in the wider region set), lets test if there is significant enrichment within
# any of our specific brain regions in the UU CRE set using a chi-squared test

# 1. Bring in the overlap counts
UU_overlap_regions <- read.delim("~/Documents/GitHub/comparativeCD/data/06_cCRE/brain_region_overlaps.txt", sep =" ")

# Calculate No Overlap
UU_overlap_regions$No_overlap <- UU_overlap_regions$Total - UU_overlap_regions$Overlap

# Extract dynamically from the dataframe so you never have to type them out
regions <- UU_overlap_regions$Region

# Extract only the two columns we need and convert to a matrix
count_matrix <- as.matrix(UU_overlap_regions[, c("Overlap", "No_overlap")])

# Assign the brain regions as the row names of the matrix
rownames(count_matrix) <- regions

# 2. Now perform the Global Chi-Square Test
global_test <- chisq.test(count_matrix)
# N.S. p = 0.44

# What if we look at number of bp overlaps, rather than element counts...

# 1. Bring in the overlap counts
#/proj/caninebrain_2025/CCD_brain_CRE_intersect/CRE_overlap/UU_elements
UU_overlaps <- read.delim("~/Documents/GitHub/comparativeCD/data/06_cCRE/brain_region_bp_overlaps.txt", sep =" ")

# Calculate No Overlap
UU_overlaps$No_overlap_bp <- UU_overlaps$Total_bp - UU_overlaps$Overlap_bp

# Extract dynamically from the dataframe so you never have to type them out
regions <- UU_overlaps$Region

# Extract only the two columns we need and convert to a matrix
count_matrix <- as.matrix(UU_overlaps[, c("Overlap_bp", "No_overlap_bp")])

# Assign the brain regions as the row names of the matrix
rownames(count_matrix) <- regions
count_matrix_dogCD  <- count_matrix
# 2. Now perform the Global Chi-Square Test
global_test <- chisq.test(count_matrix)

# 3. "One vs. Rest" Fisher's Tests with FDR Correction
# We only proceed if the global test suggests a difference exists
if(global_test$p.value < 0.05) {
  cat("\nGlobal test significant! Running Post-Hoc 'One vs Rest' tests...\n")
  
  # Create empty vectors to store our raw p-values and odds ratios
  raw_p_values <- numeric(length(regions))
  odds_ratios <- numeric(length(regions))
  
  for(i in 1:length(regions)) {
    # 'One' is the current region - pulling directly from the dataframe columns!
    one_overlap <- UU_overlaps$Overlap_bp[i]
    one_no_overlap <- UU_overlaps$No_overlap_bp[i]
    
    # 'Rest' is the sum of all OTHER regions
    rest_overlap <- sum(UU_overlaps$Overlap_bp[-i])
    rest_no_overlap <- sum(UU_overlaps$No_overlap_bp[-i])
    
    # Build the 2x2 table for this specific test
    temp_table <- matrix(c(one_overlap, one_no_overlap, 
                           rest_overlap, rest_no_overlap), 
                         nrow = 2, byrow = TRUE)
    
    # Run one-sided Fisher's test (checking for enrichment: One > Rest)
    f_test <- fisher.test(temp_table, alternative = "greater")
    
    # Store results
    raw_p_values[i] <- f_test$p.value
    odds_ratios[i] <- f_test$estimate
  }
  
  # Apply Benjamini-Hochberg FDR correction to the p-values
  fdr_adjusted_p <- p.adjust(raw_p_values, method = "BH")
  
  # Compile everything into a neat results table
  results_df <- data.frame(
    Region = regions,
    Odds_Ratio = round(odds_ratios, 2),
    Raw_P = raw_p_values,
    FDR_Adjusted_P = fdr_adjusted_p
  )
  
  # Sort by most significant to least significant
  results_df <- results_df[order(results_df$FDR_Adjusted_P), ]
  
  cat("\n--- Final Enrichment Results ---\n")
  print(results_df)
  
} else {
  cat("\nGlobal test is not significant. No single region is driving the overlap more than others.\n")
}

# Plot the results
# 1. Calculate the absolute percentage of overlapping base pairs
UU_overlaps$Percent_Overlap <- (UU_overlaps$Overlap_bp / UU_overlaps$Total_bp) * 100

# 2. Reorder the 'Region' factor so the plot sorts from highest to lowest
# This makes it instantly visually clear which region is the top hit
UU_overlaps$Region <- factor(UU_overlaps$Region, 
                             levels = UU_overlaps$Region[order(UU_overlaps$Percent_Overlap, decreasing = TRUE)])

# 3. Create the bar plot
overlap_plot <- ggplot(UU_overlaps, aes(x = Region, y = Percent_Overlap, fill = Region)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +
  
  # Add the exact percentage text on top of each bar
  geom_text(aes(label = sprintf("%.2f%%", Percent_Overlap)), 
            vjust = -0.5, size = 4.5) +
  
  # Expand the y-axis slightly so the text doesn't get cut off at the top
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  
  # Apply a clean, professional theme
  theme_minimal(base_size = 14) +
  labs(
    title = "Percentage of dogCD GWS Overlap by Brain Region",
    x = "Brain Region",
    y = "Percentage of Base Pairs Overlapping (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank() # Removes vertical grid lines for a cleaner look
  )

# 4. Display the plot
print(overlap_plot)

###SIZE:##########
# 1. Bring in the overlap bp counts

UU_overlaps <- read.delim("~/Documents/GitHub/comparativeCD/data/06_cCRE/SIZE_brain_region_bp_overlaps.txt", sep ="\t")

# Calculate No Overlap
UU_overlaps$No_overlap_bp <- UU_overlaps$Total_bp - UU_overlaps$Overlap_bp

# Extract dynamically from the dataframe so you never have to type them out
regions <- UU_overlaps$Region

# Extract only the two columns we need and convert to a matrix
count_matrix <- as.matrix(UU_overlaps[, c("Overlap_bp", "No_overlap_bp")])
#

# Assign the brain regions as the row names of the matrix
rownames(count_matrix) <- regions
count_matrix_size  <- count_matrix

# 2. Now perform the Global Chi-Square Test
global_test <- chisq.test(count_matrix)

# 3. "One vs. Rest" Fisher's Tests with FDR Correction
# We only proceed if the global test suggests a difference exists
if(global_test$p.value < 0.05) {
  cat("\nGlobal test significant! Running Post-Hoc 'One vs Rest' tests...\n")
  
  # Create empty vectors to store our raw p-values and odds ratios
  raw_p_values <- numeric(length(regions))
  odds_ratios <- numeric(length(regions))
  
  for(i in 1:length(regions)) {
    # 'One' is the current region - pulling directly from the dataframe columns!
    one_overlap <- UU_overlaps$Overlap_bp[i]
    one_no_overlap <- UU_overlaps$No_overlap_bp[i]
    
    # 'Rest' is the sum of all OTHER regions
    rest_overlap <- sum(UU_overlaps$Overlap_bp[-i])
    rest_no_overlap <- sum(UU_overlaps$No_overlap_bp[-i])
    
    # Build the 2x2 table for this specific test
    temp_table <- matrix(c(one_overlap, one_no_overlap, 
                           rest_overlap, rest_no_overlap), 
                         nrow = 2, byrow = TRUE)
    
    # Run one-sided Fisher's test (checking for enrichment: One > Rest)
    f_test <- fisher.test(temp_table, alternative = "greater")
    
    # Store results
    raw_p_values[i] <- f_test$p.value
    odds_ratios[i] <- f_test$estimate
  }
  
  # Apply Benjamini-Hochberg FDR correction to the p-values
  fdr_adjusted_p <- p.adjust(raw_p_values, method = "BH")
  
  # Compile everything into a neat results table
  results_df <- data.frame(
    Region = regions,
    Odds_Ratio = round(odds_ratios, 2),
    Raw_P = raw_p_values,
    FDR_Adjusted_P = fdr_adjusted_p
  )
  
  # Sort by most significant to least significant
  results_df <- results_df[order(results_df$FDR_Adjusted_P), ]
  
  cat("\n--- Final Enrichment Results ---\n")
  print(results_df)
  
} else {
  cat("\nGlobal test is not significant. No single region is driving the overlap more than others.\n")
}



# Plot the results
# 1. Calculate the absolute percentage of overlapping base pairs
UU_overlaps$Percent_Overlap <- (UU_overlaps$Overlap_bp / UU_overlaps$Total_bp) * 100

# 2. Reorder the 'Region' factor so the plot sorts from highest to lowest
# This makes it instantly visually clear which region is the top hit
UU_overlaps$Region <- factor(UU_overlaps$Region, 
                             levels = UU_overlaps$Region[order(UU_overlaps$Percent_Overlap, decreasing = TRUE)])

# 3. Create the bar plot
overlap_plot <- ggplot(UU_overlaps, aes(x = Region, y = Percent_Overlap, fill = Region)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +
  
  # Add the exact percentage text on top of each bar
  geom_text(aes(label = sprintf("%.2f%%", Percent_Overlap)), 
            vjust = -0.5, size = 4.5) +
  
  # Expand the y-axis slightly so the text doesn't get cut off at the top
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  
  # Apply a clean, professional theme
  theme_minimal(base_size = 14) +
  labs(
    title = "Percentage of SIZE GWAS Overlap by Brain Region",
    x = "Brain Region",
    y = "Percentage of Base Pairs Overlapping (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank() # Removes vertical grid lines for a cleaner look
  )

# 4. Display the plot
print(overlap_plot)


####################################
#### compare dogCD vs size:
#total dogCD_gws_clump100kb region encompass: 6790724
#total size_gws_clump100kb region encompass: 6377108
doggenome = 2400000000 
dogCDbp = 6790724
Size_bp = 6377108

dogCDbp/doggenome * 100
Size_bp/doggenome * 100

count_matrix_size
count_matrix_dogCD

count_bp_dogCD <- as.data.frame(count_matrix_dogCD)
sum(count_bp_dogCD$Overlap_bp)
mean(count_bp_dogCD$Overlap_bp+count_bp_dogCD$No_overlap_bp)


##testing:
dogCDbp <- 6790724
Size_bp  <- 6377108

dogCD_overlap <- as.data.frame(count_matrix_dogCD)$Overlap_bp
size_overlap  <- as.data.frame(count_matrix_size)$Overlap_bp

regions <- rownames(count_matrix_dogCD)

results <- data.frame(
  region       = regions,
  rate_dogCD   = dogCD_overlap / dogCDbp,
  rate_size    = size_overlap  / Size_bp,
  OR           = NA,
  p            = NA
)

for (i in seq_along(regions)) {
  mat <- matrix(c(dogCD_overlap[i],  dogCDbp - dogCD_overlap[i],
                  size_overlap[i],   Size_bp  - size_overlap[i]),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("dogCD","size"), c("overlap","no_overlap")))
  ft <- fisher.test(mat)
  results$OR[i] <- ft$estimate
  results$p[i]  <- ft$p.value
}

results$p_fdr <- p.adjust(results$p, method = "fdr")
print(results)


#####

dogCDbp <- 6790724
Size_bp  <- 6377108

dogCD_overlap <- as.data.frame(count_matrix_dogCD)$Overlap_bp
size_overlap  <- as.data.frame(count_matrix_size)$Overlap_bp

regions <- rownames(count_matrix_dogCD)

results <- data.frame(
  region       = regions,
  rate_dogCD   = dogCD_overlap / dogCDbp,
  rate_size    = size_overlap  / Size_bp,
  OR           = NA,
  CI_low       = NA,
  CI_high      = NA,
  p            = NA
)

for (i in seq_along(regions)) {
  mat <- matrix(c(dogCD_overlap[i],  dogCDbp - dogCD_overlap[i],
                  size_overlap[i],   Size_bp  - size_overlap[i]),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("dogCD","size"), c("overlap","no_overlap")))
  ft <- fisher.test(mat)
  results$OR[i]      <- ft$estimate
  results$CI_low[i]  <- ft$conf.int[1]
  results$CI_high[i] <- ft$conf.int[2]
  results$p[i]       <- ft$p.value
}

results$p_fdr <- p.adjust(results$p, method = "fdr")
print(results)

# Summarize the per-region ORs
mean(results$OR)   # mean OR across regions
range(results$OR)  # range


####################
# Pooled rate
###################
cat("dogCD overall rate:", sum(dogCD_overlap) / (dogCDbp * nrow(count_matrix_dogCD)), "\n")
cat("size overall rate: ", sum(size_overlap)  / (Size_bp  * nrow(count_matrix_size)),  "\n")

# Mantel-Haenszel pooled OR across regions
library(metafor)

dogCDbp <- 6790724
Size_bp  <- 6377108

dogCD_overlap <- as.data.frame(count_matrix_dogCD)$Overlap_bp
size_overlap  <- as.data.frame(count_matrix_size)$Overlap_bp

# or simply a pooled Fisher on summed counts:
mat_pooled <- matrix(c(sum(dogCD_overlap),  dogCDbp*8 - sum(dogCD_overlap),
                       sum(size_overlap),   Size_bp*8  - sum(size_overlap)),
                     nrow = 2, byrow = TRUE)
fisher.test(mat_pooled)


# Format for mantelhaen.test: needs a 2x2xK array
arr <- array(NA, dim = c(2, 2, length(regions)),
             dimnames = list(c("dogCD","size"), c("overlap","no_overlap"), regions))

for (i in seq_along(regions)) {
  arr[,,i] <- matrix(c(dogCD_overlap[i],  dogCDbp - dogCD_overlap[i],
                       size_overlap[i],   Size_bp  - size_overlap[i]),
                     nrow = 2, byrow = TRUE)
}

mantelhaen.test(arr)


###
## Two-panel figure: dogCD vs. body-size GWAS region overlap with brain cCREs
## Panel A: overlap rates (dogCD vs. body size) per brain region
## Panel B: odds ratio (95% CI), Fisher's exact test, per brain region
##
## df is the `results` object computed above (per-region Fisher's exact test,
## dogCD vs. body-size GWAS-region overlap with brain cCREs), not a hardcoded
## copy of its numbers, so this figure always reflects the live computation.

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

df <- results[, c("region", "rate_dogCD", "rate_size", "OR", "CI_low", "CI_high")]

# Order regions by OR (ascending) so the strongest enrichment sits at the top
# in both panels, and carry that order to the rates panel for alignment
region_order <- df$region[order(df$OR)]
df <- df %>%
  mutate(
    region    = factor(region, levels = region_order),
    highlight = ifelse(as.character(region) == "ACG", "ACG", "Other regions"),
    label     = sprintf("%.3f (%.3f\u2013%.3f)", OR, CI_low, CI_high)
  )

## ---- Panel A: overlap rates ----
rates_long <- df %>%
  select(region, rate_dogCD, rate_size) %>%
  pivot_longer(cols = c(rate_dogCD, rate_size), names_to = "trait", values_to = "rate") %>%
  mutate(trait = recode(trait, rate_dogCD = "dogCD", rate_size = "Body size"))

p_rates <- ggplot(rates_long, aes(x = rate, y = region, fill = trait)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  scale_fill_manual(values = c("dogCD" = "#D55E00", "Body size" = "grey65"), name = NULL) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = "Overlap rate with brain cCREs", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.margin = margin(b = -4, l = -2),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

## ---- Panel B: odds ratio forest plot ----
p_or <- ggplot(df, aes(y = region)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = CI_low, xend = CI_high, yend = region, color = highlight),
               linewidth = 1.0, lineend = "round") +
  geom_point(aes(x = OR, color = highlight), size = 3.2) +
  geom_text(aes(x = 1.43, label = label), hjust = 0, size = 2.95, family = "sans", color = "grey20") +
  scale_color_manual(values = c("ACG" = "#D55E00", "Other regions" = "#2C5F7C"), guide = "none") +
  scale_x_continuous(limits = c(1.00, 1.95), breaks = seq(1.0, 1.4, 0.1)) +
  labs(x = "Odds ratio (dogCD vs. body size)", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    plot.margin = margin(t = 23),  # nudge down to align with panel A under its legend
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

combined <- p_rates + p_or +
  plot_layout(widths = c(1, 1.25)) +
  plot_annotation(
    title = "dogCD GWAS region overlap with brain cCREs, by region",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 12, family = "sans"))
  )

ggsave("forest_plot_dogCD.pdf", combined, width = 9.5, height = 4.0, device = cairo_pdf)
ggsave("forest_plot_dogCD.png", combined, width = 9.5, height = 4.0, dpi = 300)

######
## Two-panel figure: dogCD vs. body-size GWAS region overlap with brain cCREs
## Panel A: overlap rates (dogCD vs. body size) per brain region
## Panel B: odds ratio (95% CI), Fisher's exact test, per brain region
##
## df is the `results` object computed above (per-region Fisher's exact test,
## dogCD vs. body-size GWAS-region overlap with brain cCREs), not a hardcoded
## copy of its numbers, so this figure always reflects the live computation.

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

df <- results[, c("region", "rate_dogCD", "rate_size", "OR", "CI_low", "CI_high")]

# Order regions by OR (ascending) so the strongest enrichment sits at the top
# in both panels, and carry that order to the rates panel for alignment
region_order <- df$region[order(df$OR)]
df <- df %>%
  mutate(
    region    = factor(region, levels = region_order),
    highlight = ifelse(as.character(region) == "ACG", "ACG", "Other regions"),
    label     = sprintf("%.3f (%.3f\u2013%.3f)", OR, CI_low, CI_high)
  )

## ---- Panel A: overlap rates ----
rates_long <- df %>%
  select(region, rate_dogCD, rate_size) %>%
  pivot_longer(cols = c(rate_dogCD, rate_size), names_to = "trait", values_to = "rate") %>%
  mutate(trait = recode(trait, rate_dogCD = "dogCD", rate_size = "Body size"))

p_rates <- ggplot(rates_long, aes(x = rate, y = region, fill = trait)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  scale_fill_manual(values = c("dogCD" = "#D55E00", "Body size" = "grey65"), name = NULL) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = "Overlap rate with brain cCREs", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.margin = margin(b = -4, l = -2),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

## ---- Panel B: odds ratio forest plot ----
p_or <- ggplot(df, aes(y = region)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high, color = highlight),
                orientation = "y", width = 0.45, linewidth = 0.8) +
  geom_point(aes(x = OR, color = highlight), size = 1.6) +
  geom_text(aes(x = 1.40, label = label), hjust = 0, size = 2.95, family = "sans", color = "grey20") +
  scale_color_manual(values = c("ACG" = "#D55E00", "Other regions" = "#2C5F7C"), guide = "none") +
  scale_x_continuous(limits = c(1.00, 1.72), breaks = seq(1.0, 1.4, 0.1)) +
  labs(x = "Odds ratio (dogCD vs. body size)", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    plot.margin = margin(t = 23),  # nudge down to align with panel A under its legend
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

combined <- p_rates + p_or +
  plot_layout(widths = c(1, 1.25)) +
  plot_annotation(
    title = "dogCD GWAS region overlap with brain cCREs, by region",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 12, family = "sans"))
  )

ggsave("forest_plot_dogCD.pdf", combined, width = 9.5, height = 4.0, device = cairo_pdf)
ggsave("forest_plot_dogCD.png", combined, width = 9.5, height = 4.0, dpi = 300)


###

## Three-panel figure: dogCD vs. body-size GWAS region overlap with brain cCREs
## Panel A: total cCRE size per brain region, in Mb
## Panel B: overlap rates (dogCD vs. body size) per brain region
## Panel C: odds ratio (95% CI), Fisher's exact test, per brain region
##
## df is the `results` object computed above (per-region Fisher's exact test,
## dogCD vs. body-size GWAS-region overlap with brain cCREs), not a hardcoded
## copy of its numbers, so this figure always reflects the live computation.

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

df <- results[, c("region", "rate_dogCD", "rate_size", "OR", "CI_low", "CI_high")]

# Order regions by OR (ascending) so the strongest enrichment sits at the top
# in both panels, and carry that order to the rates panel for alignment
region_order <- df$region[order(df$OR)]
df <- df %>%
  mutate(
    region    = factor(region, levels = region_order),
    highlight = ifelse(as.character(region) == "ACG", "ACG", "Other regions"),
    label     = sprintf("%.3f (%.3f\u2013%.3f)", OR, CI_low, CI_high)
  )

## ---- cCRE size per region (total bp of brain cCREs annotated to that region) ----
size_df <- data.frame(
  region   = c("ACG", "cerebellum", "frontal", "hypothalamus",
               "occipital", "striatum", "temporal", "thalamus"),
  total_bp = c(109763200, 175449000, 185010400, 152931400,
               182116400, 154493200, 147870600, 213755800)
) %>%
  mutate(
    region   = factor(region, levels = region_order),
    total_mb = total_bp / 1e6,
    mb_label = paste0(round(total_mb), " Mb")
  )

## ---- Panel A: cCRE size per region, in Mb ----
p_size <- ggplot(size_df, aes(x = total_mb, y = region)) +
  geom_col(fill = "grey75", width = 0.62) +
  geom_text(aes(x = total_mb + 8, label = mb_label), hjust = 0, size = 2.95,
            family = "sans", color = "grey20") +
  scale_x_continuous(limits = c(0, 290), expand = expansion(mult = c(0, 0))) +
  labs(x = "Brain-region cCRE size (Mb)", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    plot.margin = margin(t = 23),  # align with panel B under its legend
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

## ---- Panel B: overlap rates ----
rates_long <- df %>%
  select(region, rate_dogCD, rate_size) %>%
  pivot_longer(cols = c(rate_dogCD, rate_size), names_to = "trait", values_to = "rate") %>%
  mutate(trait = recode(trait, rate_dogCD = "dogCD", rate_size = "Body size"))

p_rates <- ggplot(rates_long, aes(x = rate, y = region, fill = trait)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  scale_fill_manual(values = c("dogCD" = "#D55E00", "Body size" = "grey65"), name = NULL) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = "Overlap rate with brain cCREs", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.margin = margin(b = -4, l = -2),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

## ---- Panel C: odds ratio forest plot ----
p_or <- ggplot(df, aes(y = region)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high, color = highlight),
                orientation = "y", width = 0.45, linewidth = 0.8) +
  geom_point(aes(x = OR, color = highlight), size = 1.6) +
  geom_text(aes(x = 1.40, label = label), hjust = 0, size = 2.95, family = "sans", color = "grey20") +
  scale_color_manual(values = c("ACG" = "#D55E00", "Other regions" = "#2C5F7C"), guide = "none") +
  scale_x_continuous(limits = c(1.00, 1.72), breaks = seq(1.0, 1.4, 0.1)) +
  labs(x = "Odds ratio (dogCD vs. body size)", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    plot.margin = margin(t = 23),  # nudge down to align with panel A under its legend
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

combined <- p_size + p_rates + p_or +
  plot_layout(widths = c(1, 1, 1.25)) +
  plot_annotation(
    title = "dogCD GWAS region overlap with brain cCREs, by region",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 12, family = "sans"))
  )

ggsave("forest_plot2_dogCD.pdf", combined, width = 12.5, height = 4.0, device = cairo_pdf)
ggsave("forest_plot2_dogCD.png", combined, width = 12.5, height = 4.0, dpi = 300)

