#!/usr/bin/env Rscript
# run_efa.R — EFA diagnostics and fit for the 14 dogCD survey items
#
# All diagnostic prints in this script (nfactors/KMO/Bartlett/det/cronbach's alpha/loadings) feed
# nothing downstream — the factor solution actually used (F1/F2/F3's item lists) is a fixed,
# pre-decided grouping, not derived programmatically here.

suppressMessages({
  library(dplyr)
  library(data.table)
  library(psych)
  library(nFactors)
  library(ltm)
})

args <- commandArgs(trailingOnly = TRUE)
response_df_path   <- args[1]
questions_csv_path <- args[2]

response_df2 <- read.table(response_df_path, header = TRUE, sep = "\t", check.names = FALSE)
response_df3 <- response_df2[-1]

sink("efa_diagnostics.log", split = TRUE)

# remove NAs to continue
# (the redundant nfactors() call on the not-yet-NA-omitted data has been removed here)
response_df4 <- na.omit(response_df3)
fa.ocd <- nfactors(response_df4, n = 10)

nrow(response_df3)
nrow(response_df4)
colnames(response_df3)

# perform parallel analysis (unseeded — see the pipeline README)
fa.ocd.ap <- parallel(subject = nrow(response_df4),
                      var = ncol(response_df4),
                      rep = 100,
                      cent = .05)

fa.ocd.nf <- nfactors(response_df4, n = 25)
fa.ocd.ev <- eigen(cor(as.matrix(response_df4)))
fa.ocd.ns <- nScree(x = fa.ocd.ev$values, aparallel = fa.ocd.ap$eigen$qevpea)

mat_cor <- cor(as.matrix(response_df4))
mat <- response_df4

#### KMO ####
fa.kmo <- KMO(r = mat_cor)
fa.kmo
fa.kmo$MSAi

#### Bartlett's test ####
fa.bart <- cortest.bartlett(mat_cor, n = nrow(mat))
fa.bart

#### determinant ####
fa.det <- det(cor(mat_cor))
fa.det

#### cronbach's alpha ####
cronbach.alpha(response_df4, CI = TRUE)

#### run factor analysis ####
nFac <- 3
fa.varimax <- fa(r = mat_cor, cor = "poly",
                 nfactors = nFac,
                 fm = "pa",
                 max.iter = 200,
                 rotate = "varimax", SMC = FALSE)

fa.varimax$loadings
fa.varimax$uniquenesses

sink()

png("parallel_analysis_plot.png", width = 7, height = 7, units = "in", res = 150)
plotParallel(fa.ocd.ap)
dev.off()

png("nscree_plot.png", width = 7, height = 7, units = "in", res = 150)
plotnScree(fa.ocd.ns)
dev.off()

# Item Uniqueness plot, labeled with survey metadata
questions <- read.csv(questions_csv_path, sep = ";")

fa.varimax.uniqueness <- data.table(question = names(fa.varimax$uniquenesses),
                                    uniqueness = fa.varimax$uniquenesses) %>%
  merge((questions %>% arrange(id) %>%
           dplyr::select(id, string, source = tags, options, survey = title) %>%
           mutate(id = as.character(id))),
        by.x = "question", by.y = "id") %>%
  arrange(as.numeric(question))

png("uniqueness_plot.png", width = 7, height = 7, units = "in", res = 150)
plot(fa.varimax.uniqueness$question, fa.varimax.uniqueness$uniqueness)
dev.off()
