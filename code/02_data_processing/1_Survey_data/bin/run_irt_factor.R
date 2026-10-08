#!/usr/bin/env Rscript
# run_irt_factor.R — per-factor IRT fit and factor scores
#
# Genericized across the three factors via arguments rather than three near-duplicate scripts.
# Real per-factor differences preserved, not normalized away:
#   - M2()'s `type` argument: F1 and F3 use "M2*", F2 uses "C2".
#   - Which fscores() method feeds the saved output: F1 and F3 use MAP (tabscores3), F2 uses EAP
#     (tabscores1) — both EAP and MAP are always computed and plotted for every factor; only the
#     one that gets joined into the saved CSV differs.
#   - F1 has one extra diagnostic scatter plot (an item's response vs. factor SE) that F2/F3
#     don't — handled by a separate process (PLOT_EXTRA_ITEM_SE), invoked only for F1, reading
#     this script's saved CSV output rather than duplicating any computation here.
#   - tracePlot/testInfoPlot/itemInfoPlot's theta_range and scale_color_brewer palette are
#     distinct per-factor values, passed through as config rather than hardcoded.
#
# response_df2 is subsetted by column from the already-recoded response_df_dogCDitems.txt
# checkpoint rather than re-derived from the raw long-format answers table: each factor's
# recoding only touches that factor's own items, so this gives the same values as re-deriving
# from raw data.
#
# Column selection by name (`item_cols`) replaces the original's per-factor hardcoded positional
# indices (e.g. F1's `x[,2:5]`, F2's `x[,2:4]`, F3's `x[,2:8]`).
#
# Plots are PNG, not PDF. The original's `_na_count_plots.pdf` held three separate plot() calls
# as three PDF pages within one file — PNG has no equivalent concept of multiple pages in one
# file, so this splits into three separate PNGs (na_count_vs_se, na_count_vs_score, score_vs_se).

suppressMessages({
  library(dplyr)
  library(data.table)
  library(mirt)
  library(lavaan)
  library(ggplot2)
  library(ggmirt)
})

args <- commandArgs(trailingOnly = TRUE)
response_df_path      <- args[1]
factor_id             <- args[2]                          # "F1", "F2", "F3"
item_cols             <- strsplit(args[3], ",")[[1]]       # e.g. c("7","153","154","155")
item_weights          <- as.numeric(strsplit(args[4], ",")[[1]])
m2_type               <- args[5]                           # "M2*" or "C2"
score_method          <- args[6]                           # "EAP" or "MAP" — which feeds the saved CSV
trace_theta_range     <- as.numeric(strsplit(args[7], ",")[[1]])
testinfo_theta_range  <- as.numeric(strsplit(args[8], ",")[[1]])
iteminfo_theta_range  <- as.numeric(strsplit(args[9], ",")[[1]])
iteminfo_palette      <- args[10]                          # "Set2" or "Set3"

response_df2_full <- read.table(response_df_path, header = TRUE, sep = "\t", check.names = FALSE)
response_df2 <- response_df2_full[, c("dog", item_cols)]
response_df3 <- response_df2[-1]

sink(paste0(factor_id, "_diagnostics.log"), split = TRUE)

head(response_df2)
fitGraded <- mirt(response_df3, 1, itemtype = "graded")
fitGraded

## get IRT estimates & plots
response_df3_named <- data.table::copy(response_df3)
setnames(response_df3_named, old = item_cols, new = paste0("q", item_cols))

model <- paste0("f1 =~ ", paste0("q", item_cols, collapse = "+"))
fitCTT <- cfa(model, data = response_df3_named)

# IRT solution
summary(fitGraded)

# CTT solution
standardizedsolution(fitCTT) %>%
  filter(op == "=~") %>% dplyr::select(rhs, F1 = est.std)

# IRT parameters
params <- coef(fitGraded, IRTpars = TRUE, simplify = TRUE)
round(params$items, 2)

# M2() reliably errors here ("too few degrees of freedom") for all three factors — the original
# script's own comment on this exact call documents hitting the same error every time. Run
# interactively, that error just prints and the author moved on to the next line; a batch
# Rscript run would otherwise halt the whole script on it. tryCatch preserves the former.
tryCatch(
  print(M2(fitGraded, type = m2_type, calcNULL = FALSE, na.rm = TRUE)),
  error = function(e) cat("Error:", conditionMessage(e), "\n")
)

itemfit(fitGraded, na.rm = TRUE)

sink()

n_items <- length(item_cols)

png(paste0(factor_id, "_trace_plot.png"), width = 7, height = 7, units = "in", res = 150)
print(tracePlot(fitGraded, theta_range = trace_theta_range) + labs(color = "Answer Options"))
dev.off()

png(paste0(factor_id, "_iteminfo_facet_plot.png"), width = 7, height = 7, units = "in", res = 150)
itemInfoPlot(fitGraded, seq_len(n_items), facet = TRUE)
dev.off()

png(paste0(factor_id, "_testinfo_plot.png"), width = 7, height = 7, units = "in", res = 150)
testInfoPlot(fitGraded, theta_range = testinfo_theta_range, adj_factor = .5)
dev.off()

png(paste0(factor_id, "_iteminfo_plot.png"), width = 7, height = 7, units = "in", res = 150)
print(itemInfoPlot(fitGraded, seq_len(n_items), theta_range = iteminfo_theta_range, facet = FALSE, legend = TRUE) +
        scale_color_brewer(palette = iteminfo_palette))
dev.off()

sink(paste0(factor_id, "_diagnostics.log"), append = TRUE, split = TRUE)

# create the theta values for each factor
summary(fitGraded)

sink()

tabscores1 <- as.data.frame(fscores(fitGraded, method = 'EAP', full.scores.SE = TRUE, item_weights = item_weights))
tabscores3 <- as.data.frame(fscores(fitGraded, method = 'MAP', full.scores.SE = TRUE, item_weights = item_weights))

png(paste0(factor_id, "_eap_map_plot.png"), width = 7, height = 7, units = "in", res = 150)
plot(tabscores1$F1, tabscores1$SE_F1, xlab = "Factor score", ylab = "Standard error",
     main = paste0(factor_id, "_CCD3F"))
points(tabscores3$F1, tabscores3$SE_F1, col = "blue", pch = 4, cex = .8)
legend("topright", legend = c("EAP", "MAP"), pch = 19, bty = "n", col = c("green", "blue"))
dev.off()

# combine with dog ID — using EAP for F2, MAP for F1/F3 (see header note)
chosen_scores <- if (identical(score_method, "EAP")) tabscores1 else tabscores3
chosen_scores$number <- rownames(chosen_scores)
response_df2$number <- rownames(response_df2)

sink(paste0(factor_id, "_diagnostics.log"), append = TRUE, split = TRUE)

nrow(chosen_scores)
nrow(response_df2)
result <- left_join(response_df2, chosen_scores, by = "number")
head(result)
nrow(result)

sink()

out <- data.frame(dog = result$dog, score = result$F1, score_SE = result$SE_F1)
for (col in item_cols) {
  out[[paste0("item", col)]] <- result[[col]]
}
names(out)[names(out) == "score"]    <- factor_id
names(out)[names(out) == "score_SE"] <- paste0(factor_id, "_SE")

sink(paste0(factor_id, "_diagnostics.log"), append = TRUE, split = TRUE)
head(out)
nrow(out)
sink()

x <- result
x$na_count <- apply(x[, item_cols], 1, function(row) sum(is.na(row)))

sink(paste0(factor_id, "_diagnostics.log"), append = TRUE, split = TRUE)
length(x$na_count)
length(x$F1)
sink()

png(paste0(factor_id, "_na_count_vs_se_plot.png"), width = 7, height = 7, units = "in", res = 150)
plot(x$SE_F1, x$na_count, xlab = "SE", pch = 19, ylab = "NA count", main = paste0(factor_id, "_CCD3F"))
dev.off()

png(paste0(factor_id, "_na_count_vs_score_plot.png"), width = 7, height = 7, units = "in", res = 150)
plot(x$F1, x$na_count, xlab = "Factor score", pch = 19, ylab = "NA count", main = paste0(factor_id, "_CCD3F"))
dev.off()

png(paste0(factor_id, "_score_vs_se_plot.png"), width = 7, height = 7, units = "in", res = 150)
plot(x$F1, x$SE_F1, xlab = "Factor score", pch = 19, ylab = "SE", main = paste0(factor_id, "_CCD3F"))
dev.off()

sink(paste0(factor_id, "_diagnostics.log"), append = TRUE, split = TRUE)
head(out)
nrow(out)
sink()

write.table(out, file = paste0(factor_id, "_CCD3F.txt"), row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)
