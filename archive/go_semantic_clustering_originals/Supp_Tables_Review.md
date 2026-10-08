# Supplementary Tables Review — Canine GWAS

---

## Index sheet

**Issues:**
- Row 3 (Table 3 description) contains a private working note: *"NOTE! This incorporates Sup Table 6 data - check to see if you want to keep that still, it may be that you want a version of St6 for the cCRE"* — must be removed before submission.
- Tables 4 and 5 are missing from the index (numbering jumps from 3 to 6). Either the missing tables need entries or the numbering needs to be reconciled.
- Multiple blank rows scattered throughout — clean up.
- Sheet21 is empty and should be deleted.

---

## Table 1 — Factor Analysis Loading

**Verdict: trim significantly.**

The core data (Factor 1/2/3 loadings, Item Number, Item Description) is fine and should be kept. However, columns 6–12 (tags, style, format, options, survey, title, intro) appear to be raw metadata exported directly from the survey software. These are internal database fields, not interpretable by a reader, and should be removed entirely. In particular, the "intro" column contains paragraph-length survey section blurbs that make each row enormous and add no scientific information.

**Keep:** Factor 1, Factor 2, Factor 3, Item Number, Item Description.  
**Remove:** tags, style, format, options, survey, title, intro.

---

## Table 2 — Demographics of GWAS datasets

**Verdict: fine as-is.** Clean, readable, appropriately sized. No major changes needed. Minor option: group rows by phenotype type (Continuous first, then Categorical) with a divider or subtle shading to aid scanning.

---

## Table 3 — Annotated GWAS SNPs

**Verdict: major cleanup needed.**

This is 89 columns × 843 rows. Several problems:

**1. ~28 trailing empty columns.** Delete them all.

**2. Duplicated tissue-type column headers (two sets).** The columns "Cerebellum, Cerebrum, Colon, Kidney, Liver, Lung, Mammary Gland, Ovary, Pancreas, Spleen, Stomach" appear twice consecutively, and "Anterior Cingulate Cortex, Cerebellum, Frontal Lobe, Hypothalamus, Occipital Cortex, Striatum, Temporal Cortex, Thalamus" also appears twice. It is not clear what differentiates the two sets (e.g., are they two different assays, annotation tools, or species?). Add a two-row grouped header (e.g., a parent row labeling each block with what it represents), or rename columns to make the distinction explicit (e.g., "Dog_Cerebellum" vs "Human_Cerebellum", or "H3K27ac_Cerebellum" vs "ATAC_Cerebellum").

**3. Confusing column header encoding.** "SuSiE Region [NO_PRIORS = 1; NO_PRIORS_clump+-100kb = 2]" embeds a legend inside a column name. Move this explanation to a table footnote and use a simpler column header like "SuSiE Region Type."

**4. Column "Genes to NetcoLoc From Region"** — "NetcoLoc" should probably be spelled consistently with how it appears in the manuscript.

---

## Table 6 — GWAS associated regions, credible sets, and implicated genes

**Verdict: moderate cleanup needed.**

**1. Duplicated column names.** "POS" appears in columns 6 and 8 (different content — one is the credible set range as a string "chr1:75142353..75150007", the other is just chromosome label). Rename them clearly. Similarly "GWAS" appears in columns 4 and 11 with unclear distinction.

**2. Compound data in column 6 ("POS").** Values like "chr2:57876761..57934887" embed three pieces of information (chr, start, stop) as a string. Since start and stop already exist as separate columns (10, 11), this column is redundant — remove it, or split it.

**3. Date-stamps in column headers.** Column names ending in "_251121" (presumably a date: 25 Nov 2021) are internal tracking artifacts. Remove before publication.

**4. Long notes embedded in headers.** "NO_PRIORS_clump: start stop (Nsnps CS)_251121" is a column header that requires decoding. Rename to something interpretable and move explanations to footnotes.

**5. Column "Region_sorted by chr_pos_sign"** is just a sort-order index and the name encodes the sort logic. If it's just a region ID, name it "Region ID." If it's not needed, remove it.

---

## Table 7 — Genes taken forward for network analysis

**Verdict: mostly fine, minor fixes.**

- Typo in title: "forawrd" → "forward."
- "Species" column is "Dog" for every row — redundant. Remove it or repurpose it to include human orthologs where relevant.
- "dogCD significance" is always either "significant" or (presumably) "not significant" — consider replacing with a boolean (Yes/No) or grouping rows by significance.
- "dogCD GWAS region" contains numeric region IDs with no legend explaining what region 1, 2, 4, 5, etc. correspond to. Add a footnote cross-referencing Table 6 region IDs.

---

## Table 8 — Network gene lists (OCD/SCZ/MD overlap)

**Verdict: structural rethink recommended.**

The table places six independent gene lists of different lengths side-by-side, with NaN padding wherever a gene isn't in a given list. This makes it look like row relationships exist across columns when they don't — a reader might think gene on row 5 of column A corresponds to gene on row 5 of column B.

**Recommended restructure:** Convert to a single tidy table with one row per unique gene and binary (Yes/No or 1/0) columns indicating membership in each gene set, plus retain the "Shared/not-shared" and "Illustration" (pathway) annotation columns. This is far more queryable and unambiguous.

If the side-by-side list format is kept for visual reasons, add a prominent note clarifying that rows are not linked across columns.

---

## Table 9 — Colocalization network gene lists

**Verdict: same structural problem as Table 8.**

Six independent gene lists (N = 184, 165, 462, 452, 387, 382) placed in side-by-side columns with NaN padding. The longest list has 462 rows, creating a very tall table where most cells are empty. Same recommendation as Table 8: restructure as a presence/absence matrix. The column headers make the comparisons clear (OCD_dogCD, SCZ_dogCD, MD_dogCD × network vs. systemsmap), so a matrix format would be natural here.

---

## Table 10 — Heritability estimates

**Verdict: fine, minor improvements possible.**

The data is organized well (rows are items/factors, sorted by some criterion). Improvements:

- Column headers are very long (e.g., "GREML_LDS_constrained_heritability"). Consider a two-row grouped header: parent row for "GREML constrained / GREML LDS constrained / GREML unconstrained / GREML LDS unconstrained," child row for "h² / SE" under each. This halves the number of visible column headers and makes the table scannable.
- NaN cells (where a method wasn't applicable) should use "—" instead for a cleaner published appearance.

---

## Table 11 — Sensitivity analysis (colocalization of dogCD and OCD network genes)

**Verdict: fine, one formatting issue.**

Numbers use European decimal convention (commas as decimal separators: "346,72", "7,94E-05"). For an English-language journal, these should be converted to periods ("346.72", "7.94E-05"). Check whether this was introduced by a locale setting in the analysis software.

---

## Table 12 — Control trait correlations

**Verdict: two separate tables merged into one — split them.**

The left half (columns 0–6: OCD sub-cohort, Trait, rg, SE, Z, P, P_fdr) is a genetic correlation table. The right half (columns 7–9: Trait, SNP h2, Gene sets from) is a completely different table listing trait heritabilities and reference citations. These are only loosely related and should be separate supplementary tables.

Additionally:
- The "Gene sets from" column contains full-length reference strings (entire citation text). These belong in footnotes, not table cells.
- European decimal format (same issue as Table 11): "0,36", "4,65E-40" → "0.36", "4.65E-40."
- Note: "Reumathoid Arthritis" should be "Rheumatoid Arthritis."

---

## Table 13 — cCRE overlaps

**Verdict: two sub-tables in one sheet — split or clearly delineate.**

The sheet contains two separate analyses:
1. "EPIC dog cCREs: brain vs. other tissues" (2 trait rows) — a simple comparison table.
2. "Dog brain cCREs: overlap with GWAS regions by brain region" (8 brain region rows) — a more detailed table.

These are separated only by a blank row, making it easy to miss the second table entirely. Either place them in separate sheets or use clear visual separation (e.g., a full-width bold section header row with background color, and a note in the Index sheet).

The second sub-table also has three columns at the end labeled "Enrichment OR," "CI (low)," "CI (high)" with no indication of what they represent relative to the preceding dogCD and dog size columns. Clarify what this enrichment is (enrichment of dog brain cCREs relative to a null? one trait relative to the other?).

---

## Summary of priority fixes

| Priority | Issue |
|---|---|
| Must fix | Remove working note from Index |
| Must fix | Table 11, 12: European decimal separators (commas → periods) |
| Must fix | Table 12: "Rheumatoid" spelling |
| Must fix | Table 7: "forawrd" → "forward" in title |
| Must fix | Table 1: Remove survey software metadata columns (tags, style, format, options, survey, title, intro) |
| Should fix | Table 3: Explain/rename duplicated tissue column blocks; delete ~28 empty trailing columns |
| Should fix | Table 6: Remove date-stamps from headers, rename duplicate columns, clarify compound-data column |
| Should fix | Table 12: Split into two tables; move citations to footnotes |
| Should fix | Table 13: Split two sub-tables; clarify "Enrichment OR" columns |
| Consider | Tables 8 & 9: Restructure parallel gene lists as presence/absence matrix |
| Consider | Table 10: Use grouped two-row header for heritability columns |
| Delete | Sheet21 (empty) |
