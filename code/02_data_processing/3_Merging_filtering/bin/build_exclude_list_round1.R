#!/usr/bin/env Rscript
# build_exclude_list_round1.R — derive the first-round sample exclusion list (low depth + known
# duplicates) from frozen, already-decided inputs.
#
# This is a mechanical derivation, not a human judgment call itself — the judgment already
# happened upstream (the 0.3x depth threshold, and which samples are "duplicates to keep" in the
# frozen file below). The result can be diffed against the already-committed
# data/02_data_processing/list_exclude.txt.
#
# The `fam` join key is `IID` — despite this file's `IID` column looking like a short "dogID"
# rather than a raw sequencing sample ID, joining by `IID` as written is correct here.

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
fam_qc3modi_path <- args[1]
depth_path       <- args[2]
duplicates_path  <- args[3]
fam_qc4_path     <- args[4]

fam <- read.csv(fam_qc3modi_path, sep = "\t", header = TRUE)

DA_all <- read.csv(depth_path, sep = "\t", header = FALSE, col.names = c("sampleID", "meanDepthALL"))
DA_all$IID <- as.character(DA_all$sampleID)
DA_merged <- left_join(fam, DA_all, by = "IID")

sink("build_exclude_list_round1.log", split = TRUE)

nrow(subset(DA_merged, is.na(DA_merged$meanDepthALL)))

sink()

DAdup <- read.csv(duplicates_path, sep = "\t", header = FALSE,
                   col.names = c("Sample_run_id", "depth", "sampleID", "keep_or_exclude", "IID", "dogID", "extra_column"))
DA_merged2 <- left_join(DA_merged, DAdup, by = "IID")

sink("build_exclude_list_round1.log", append = TRUE, split = TRUE)

subset(DA_merged2, !is.na(DA_merged2$keep_or_exclude))
nrow(DA_merged2)

sink()

DA_merged3 <- DA_merged2 %>%
  mutate(DEPTH_filter = ifelse(meanDepthALL < 0.3, 0, 1))

sink("build_exclude_list_round1.log", append = TRUE, split = TRUE)
head(DA_merged3)
table(DA_merged3$DEPTH_filter, exclude = NULL)
sink()

DA_merged4 <- DA_merged3 %>%
  mutate(DEPTH_filter2 = ifelse(is.na(DEPTH_filter), 1, DEPTH_filter))

sink("build_exclude_list_round1.log", append = TRUE, split = TRUE)
head(DA_merged4)
table(DA_merged4$DEPTH_filter2, exclude = NULL)
sink()

DA_merged5 <- subset(DA_merged4, DA_merged4$DEPTH_filter2 == "1")

sink("build_exclude_list_round1.log", append = TRUE, split = TRUE)
nrow(DA_merged5)
sink()

DA_merged6 <- DA_merged5 %>%
  mutate(keep_or_exclude2 = ifelse(is.na(keep_or_exclude), 1, keep_or_exclude))

sink("build_exclude_list_round1.log", append = TRUE, split = TRUE)
table(DA_merged6$keep_or_exclude2, exclude = NULL)
sink()

DA_merged7 <- subset(DA_merged6, !DA_merged6$keep_or_exclude2 == "exclude")

sink("build_exclude_list_round1.log", append = TRUE, split = TRUE)
nrow(DA_merged7)
table(DA_merged7$keep_or_exclude2, exclude = NULL)
summary(DA_merged7$meanDepthALL, exclude = NULL)
head(DA_merged7)
sink()

DA_merged7$qc_pass <- 1

famQC4 <- read.csv(fam_qc4_path, sep = " ", header = FALSE, col.names = c("FID", "IID", "F", "M", "sex", "phe"))

sink("build_exclude_list_round1.log", append = TRUE, split = TRUE)
head(famQC4)
sink()

QC4updated <- left_join(famQC4, DA_merged7, by = "FID")

sink("build_exclude_list_round1.log", append = TRUE, split = TRUE)
head(QC4updated)
table(QC4updated$qc_pass, exclude = NULL)
sink()

QC4updated2 <- QC4updated %>%
  mutate(qc_pass2 = ifelse(is.na(qc_pass), 2, qc_pass))

famQC4_exclude <- subset(QC4updated2, QC4updated2$qc_pass2 == "2")

sink("build_exclude_list_round1.log", append = TRUE, split = TRUE)
nrow(famQC4_exclude)
famQC4_exclude$FID
sink()

list_exclude <- data.frame(FID = famQC4_exclude$FID, IID = famQC4_exclude$IID.x)
write.table(list_exclude, file = "list_exclude.txt", row.names = FALSE, sep = "\t", col.names = FALSE, quote = FALSE)
