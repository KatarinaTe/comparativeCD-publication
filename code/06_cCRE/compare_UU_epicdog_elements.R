# Fig. 4d: element-count comparison of UU vs. EpicDog cCREs (Promoters/Enhancers/Repressed),
# by Shared / Unique to EPIC / Unique to UU.
#
# Based on a collaborator draft. Two fixes applied before this was run for real:
#
# 1. read_bed() crashes on these files as-is: cerebellum_13_dense.bed/cerebrum_13_dense.bed
#    (and the UU input) carry a UCSC "track name=..." header line on row 1, which
#    valr::read_bed()'s internal field-sniffing misreads, producing
#    "expected 3 required names, missing: start and end". read_dense_bed() below strips
#    that line first. Confirmed against a real sample of the fetched EpicDog format.
#
# 2. Unique_Epic/Unique_UU previously used bed_subtract() with its default any = FALSE,
#    which returns leftover coordinate FRAGMENTS after removing the overlapping portion of
#    each element rather than a whole-element count -- an element only partially overlapped
#    by the other dataset gets split into 2+ fragments and over-counted. Demonstrated with a
#    synthetic 2-element test (one partially overlapped, one not): the old method returned
#    Unique(3) + Shared(1) = 4 for a true total of 2. Since UU and EpicDog are two
#    independently-trained ChromHMM models, state boundaries essentially never land on the
#    same base pair, so partial (not full-containment) overlap is the norm here, not an edge
#    case -- this wasn't a rounding-level effect. Fixed with bed_subtract(..., any = TRUE)
#    (per the manuscript author's own suggestion, 2026-09-21) -- valr's built-in equivalent
#    of bedtools subtractBed's -A flag, dropping a whole x element the moment it overlaps
#    anything in y, in one call. Confirmed to give numerically identical Unique_Epic/
#    Unique_UU counts to an earlier two-call bed_intersect()-based workaround this repo used
#    before finding any = TRUE. Matches what the figure's own legend and axis ("Element
#    count (K)") already assume: a partition of discrete elements, not a base-pair
#    computation.
#
# State groupings (Promoters/Enhancers/Repressed for both EPIC's 13-state and UU's
# 9-state legends) are unchanged from the collaborator's draft -- left as their call, not
# re-derived here.
#
# Provenance of all_tissues_filtered.bed confirmed directly by the manuscript author
# (2026-09-21) and reproduced byte-for-byte by build_all_tissues_filtered.sh from the
# already-fetched per-region UU BEDs -- see that script's own header for detail.
# cerebellum_13_dense.bed/cerebrum_13_dense.bed are fetchable via
# analyses/06_cCRE/fetch-epicdog-chromatin-states-bed.sh.
#
# The two fixes above (read_bed() crash, bed_subtract() over-counting) were reviewed
# against the real data with the manuscript author, who approved re-plotting Fig. 4d with
# the corrected (boolean any-overlap) counts -- see REPRODUCIBILITY_AUDIT.md's 2026-09-21
# update for the before/after numbers and the reasoning.

library(tidyverse)
library(valr)

options(scipen = 10000)

# --- 1. Read inputs -----------------------------------------------------------------

read_dense_bed <- function(path) {
  lines <- readr::read_lines(path)
  if (grepl("^track", lines[1])) lines <- lines[-1]
  tmp <- tempfile(fileext = ".bed")
  readr::write_lines(lines, tmp)
  x <- read_bed(tmp)
  unlink(tmp)
  x
}

# Epic files (column 4 = state, read as character; filters below compare against numeric
# state IDs and rely on R's normal numeric<->character coercion in %in%, confirmed to work)
epic_cerebellum <- read_dense_bed("cerebellum_13_dense.bed")
epic_cerebrum   <- read_dense_bed("cerebrum_13_dense.bed")

# UU data
# TODO: confirm provenance of all_tissues_filtered.bed before running this for real.
uu_all <- read_dense_bed("all_tissues_filtered.bed")

# --- 2. Helper: filter by state, pool, and classify shared/unique whole elements -----

process_and_compare <- function(epic_cbl, epic_crb, uu_data, epic_states, uu_states, label) {

  epic_combined <- bind_rows(
    epic_cbl %>% filter(name %in% epic_states),
    epic_crb %>% filter(name %in% epic_states)
  ) %>%
    bed_merge()

  uu_merged <- uu_data %>%
    filter(name %in% uu_states) %>%
    bed_merge()

  # "Shared" reported from EPIC's own element boundaries, same convention as the
  # collaborator's original draft.
  shared_count <- bed_intersect(epic_combined, uu_merged) %>%
    distinct(chrom, start.x, end.x) %>%
    nrow()

  # bed_subtract(..., any = TRUE) drops a whole x element the moment it overlaps anything
  # in y -- the correct whole-element "unique" count directly, one call each side.
  unique_epic <- bed_subtract(epic_combined, uu_merged, any = TRUE) %>% nrow()
  unique_uu   <- bed_subtract(uu_merged, epic_combined, any = TRUE) %>% nrow()

  data.frame(
    Feature = label,
    Unique_Epic = unique_epic,
    Unique_UU = unique_uu,
    Shared = shared_count
  )
}

# --- 3. Run for each category (state groupings unchanged from the collaborator draft) --

# Promoters: Epic (6-9) vs UU (1,4,5,6)
promoters_res <- process_and_compare(epic_cerebellum, epic_cerebrum, uu_all,
                                      c(6,7,8,9), c(1,4,5,6), "Promoters")

# Enhancers: Epic (11-13) vs UU (2,8)
enhancers_res <- process_and_compare(epic_cerebellum, epic_cerebrum, uu_all,
                                      c(11,12,13), c(2,8), "Enhancers")

# Repressed: Epic (1-3) vs UU (9)
repressed_res <- process_and_compare(epic_cerebellum, epic_cerebrum, uu_all,
                                      c(1,2,3), c(9), "Repressed")

# --- 4. Final summary table and plot -------------------------------------------------

final_summary <- bind_rows(promoters_res, enhancers_res, repressed_res)
print(final_summary)

final_summary_long <- final_summary %>%
  pivot_longer(
    cols = c(Unique_Epic, Unique_UU, Shared),
    names_to = "Category",
    values_to = "Count"
  )

final_summary_long$Category <- gsub("_", " ", final_summary_long$Category)
final_summary_long$Feature <- factor(final_summary_long$Feature, levels = c("Promoters", "Enhancers", "Repressed"))

write.table(final_summary_long, file = "EPIC_comparison.txt", quote = FALSE, row.names = FALSE, sep = "\t")

UU_epic_comparison_barplot <- final_summary_long %>%
  ggplot(aes(y = Count / 1000, x = Feature, fill = Category)) +
  geom_col() +
  ylab("Element count (K)") +
  xlab("Element type") +
  scale_fill_manual(values = c("purple3", "orange3", "cyan3")) +
  theme_minimal() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

ggsave("Fig4d_corrected.pdf", UU_epic_comparison_barplot, width = 4.5, height = 6)
