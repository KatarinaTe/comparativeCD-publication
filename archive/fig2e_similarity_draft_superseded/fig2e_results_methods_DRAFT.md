# Fig. 2e — draft Results + Methods text (2026-10-01)

**Status: DRAFT, not verified against an original source script.** KatarinaTe asked to reconstruct
this analysis from the manuscript Results text alone, since the original script/notebook behind
the published Fig. 2e hasn't been located yet. This text assumes that reconstruction is correct
for now — review before using, especially the OCD–Depression finding, which is a real change from
the currently published number (see "What changed" below).

Full computation: `code/07_additional_plots_and_analyses/Figure2_new_panels/` session work,
script and outputs saved to
`data/results/05_network/pipeline_output/fig2/fig2e_similarity_heatmap_DRAFT.{R,png,pdf}` and
`fig2e_similarity_results_DRAFT.csv` (local-only, not committed — `data/results/` is gitignored).
**Not yet wired into the Nextflow pipeline** — planned for after KatarinaTe's review.

---

## Results text (draft)

> Measuring semantic similarity across all enriched terms confirms this (Fig. 2e): the OCD/dogCD
> network is highly similar to both the depression and schizophrenia networks (0.77 and 0.66 of
> the maximum possible similarity above chance; permutation P = 0.001 for each), more similar than
> depression and schizophrenia are to each other (0.50; P = 0.001). The human-only OCD network, by
> contrast, shows no greater similarity than chance to depression (0.00; P = 0.49) and is less
> similar than random term sets to schizophrenia (−0.13) or to the OCD/dogCD network (−0.14).

## Methods text (draft)

> *Semantic similarity across enriched GO terms (Fig. 2e).* For each of the four gene-set
> hierarchies compared in Fig. 2 (OCD, human-only; OCD/dogCD, cross-species; Depression,
> human-only; Schizophrenia, human-only), we defined an enriched-term set as all Gene Ontology
> Biological Process terms reaching nominal significance (−log₁₀P ≥ 5) in that hierarchy's
> GO-term enrichment analysis, giving N = 8 (OCD), 75 (OCD/dogCD), 129 (Depression) and 116
> (Schizophrenia) terms. For every pair of hierarchies, we computed the semantic similarity
> between their term sets using the Wang method (GOSemSim v2.26.1 in R), combined across all term
> pairs via the best-match average (BMA). To express each pairwise similarity as a fraction of the
> maximum possible similarity achievable above chance, we built a null distribution from 1,000
> permutations per pair, each permutation drawing two random term sets of the same sizes from the
> full set of GO Biological Process terms annotated in the org.Hs.eg.db genome-wide annotation
> (n = 12,588 terms) and recomputing the Wang/BMA similarity. The chance-corrected similarity was
> calculated as (observed similarity − mean null similarity) / (1 − mean null similarity); a
> one-sided permutation P value was calculated as the proportion of null-distribution values
> greater than or equal to the observed similarity.

---

## What's solid vs. reconstructed

**Solid (deterministic, fully confirmed against the published panel):**
- Term-set definitions and sizes (N = 8/75/129/116) match the panel's own labels exactly.
- The `(N)` shown under each score in the published panel is the literal GO-term overlap count
  between each pair's term sets — verified to match exactly for all 6 pairs (0, 7, 0, 66, 46, 50).

**Reconstructed (not verified):**
- The GOSemSim Wang/BMA similarity measure itself — correctly reproduces the *ranking* of all 6
  published pairs.
- The chance-correction procedure (permutation background = full GO:BP ontology, 1,000
  permutations, score = (raw − chance)/(1 − chance)). This was chosen after an earlier attempt
  (permuting within the paper's own already-enriched term pool) gave worse and in one case
  wrong-signed results. The 3 strongest pairs land on empirical P = 0.001 exactly, matching the
  manuscript's stated P-values, which is reassuring but not proof this is the real method.

## What changed from the currently published numbers

| pair | published | this draft |
|---|---|---|
| OCD/dogCD–Depression | 0.79* | 0.77* |
| OCD/dogCD–Schizophrenia | 0.67* | 0.66* |
| Depression–Schizophrenia | 0.49* | 0.50* |
| **OCD–Depression** | **0.34\*, P = 0.001** | **0.00, P = 0.49 (not significant)** |
| OCD–OCD/dogCD | −0.23 | −0.14 |
| OCD–Schizophrenia | −0.26 | −0.13 |

The OCD–Depression cell is the one substantive change: the published claim that human-only OCD is
"modestly similar" to depression doesn't hold up under this reconstruction — OCD's 8 significant
terms are all generic transcription-regulation GO terms (not neuro/synaptic, unlike the other
three groups), which may explain the discrepancy. Worth a second pair of eyes, or the real source
script, before this replaces the published text.

---

## Exploratory extension: adding dogCD as a 5th group (2026-10-02)

KatarinaTe asked to see what the panel looks like with the single-species dogCD-only network
added as a 5th group (not part of the original published Fig. 2e, which only has 4 groups).
**Pausing here — she's discussing with a collaborator before deciding whether to pursue this
further.**

dogCD's term set was derived the same way as the other 4 groups (top-level community + −log₁₀P ≥ 5)
from `data/05_network/CCD_ONLY_hierachy_full_GO_enrichment.tsv` (top community C257, 146 terms
tested, 8 significant). The 6 original pairs are unchanged; only the 4 new dogCD pairs were
computed (same 1,000-permutation methodology):

| pair | score | overlap (N) | P |
|---|---|---|---|
| dogCD–OCD/dogCD | 0.27* | 5 | 0.001 |
| dogCD–Depression | 0.22* | 5 | 0.001 |
| dogCD–Schizophrenia | 0.17* | 5 | 0.001 |
| dogCD–OCD (human-only) | −0.03 | 0 | not significant |

Reading: dogCD alone is weakly-but-significantly similar to all three "real" neuropsychiatric
networks (Depression, Schizophrenia, and its own cross-species descendant OCD/dogCD), well below
OCD/dogCD's much stronger 0.77/0.66 — consistent with a story where dogCD alone carries some
shared signal, and combining it with human OCD data amplifies that considerably. dogCD and
human-only OCD remain unrelated, consistent with human-only OCD being the outlier throughout this
panel.

Outputs: `data/results/05_network/pipeline_output/fig2/fig2e_similarity_heatmap_DRAFT_with_dogCD.{R,png,pdf}`
and `fig2e_similarity_results_DRAFT_with_dogCD.csv` — a separate set of files, the original
4-group outputs (above) were left untouched. Same caveats apply as the rest of this document
(reconstructed normalization, not independently verified).

---

## Real legend received from collaborator (2026-10-03) — confirms the method, flags one open question

KatarinaTe's collaborator independently tried the same reconstruction and sent back the actual
figure legend intended for this panel. It confirms nearly everything in this reconstruction was
already correct:

- Significance cutoff P < 1×10⁻⁵ (−log₁₀P ≥ 5) — matches.
- "Similarity is the average of the two directions, each the mean over one network's terms of the
  highest Wang semantic similarity to any term in the other network" — this is exactly GOSemSim's
  own BMA (best-match average) definition, already what was used.
- Score formula (observed − expected)/(1 − expected), expected from 1,000 random-term-set pairs —
  matches exactly.
- Parenthetical cell counts = shared significant-term counts — already independently verified.
- Network definitions (OCD/dogCD cross-species; dogCD single-species dog-only; OCD/Depression/
  Schizophrenia human-only) — matches, and confirms the **real published figure includes dogCD as
  a 5th group**, same as the exploratory extension above.

**One thing this reconstruction was missing**: the legend specifies asterisks come from
**Benjamini-Hochberg-corrected** empirical P values across all pairs, not raw per-pair empirical P.
Fixed in the v3 rebuild below.

**Still not specified by the legend, and the one known point of disagreement**: what background
pool the 1,000 random GO term sets are drawn from. KatarinaTe's collaborator ran their own
implementation and got numbers matching this reconstruction closely on every pair **except
OCD–Depression** — and the discrepancy was stable across two of their runs (before and after adding
dogCD), so it isn't permutation noise; it's a genuine implementation difference. Likely candidate:
this reconstruction samples from the full genome-wide GO:BP annotation (org.Hs.eg.db, n = 12,588
terms), and the collaborator's implementation may use a different background (e.g. the enrichment
analysis's own tested-term universe rather than the whole ontology) — unconfirmed, pending their
exact code or their raw/expected values for that one cell.

### Legend (collaborator's wording, used as-is — already generically covers all 5 groups)

> Semantic similarity between the GO Biological Process terms significantly enriched (P < 1 × 10⁻⁵)
> in each pair of networks, relative to chance. Similarity is the average of the two directions,
> each the mean over one network's terms of the highest Wang semantic similarity to any term in the
> other network, so that both networks count equally regardless of size. Values give the fraction
> of the possible above-chance similarity achieved, (observed − expected) / (1 − expected), where
> expected is the mean for 1,000 pairs of random GO term sets of the same sizes (0, as similar as
> random terms; 1, identical sets; negative, less similar than random). Axis labels give each
> network's total number of significant GO terms; numbers in parentheses within cells give the
> number of significant GO terms shared by the two networks. *, more similar than random
> (Benjamini-Hochberg-corrected empirical P < 0.05). OCD/dogCD, colocalized cross-species network
> for OCD and dogCD; dogCD, single-species, dog-only network; OCD, depression and schizophrenia,
> human-only networks.

**Short legend:**

> Semantic similarity of enriched GO terms between networks, relative to chance. Values are
> (observed − expected)/(1 − expected), where expected is the mean for 1,000 random GO term sets of
> matching size (0, chance; 1, identical). Axis labels give significant terms per network; numbers
> in parentheses within cells give shared terms. *P < 0.05, Benjamini-Hochberg corrected.

### Methods text (v3, rebuilt against this legend)

> *Semantic similarity across enriched GO terms (Fig. 2e).* For each of five gene-set hierarchies
> (OCD, human-only; OCD/dogCD, the colocalized cross-species network for OCD and dogCD; dogCD,
> single-species, dog-only; Depression, human-only; Schizophrenia, human-only), we defined an
> enriched-term set as all Gene Ontology Biological Process terms reaching significance at
> P < 1 × 10⁻⁵ in that hierarchy's GO-term enrichment analysis, giving N = 8 (OCD), 75 (OCD/dogCD),
> 8 (dogCD), 129 (Depression) and 116 (Schizophrenia) terms. For every pair of hierarchies, we
> computed semantic similarity between their term sets using the Wang method (GOSemSim v2.26.1,
> org.Hs.eg.db v3.17.0, R), combined across both directions via the best-match average (BMA): the
> mean, over one network's terms, of each term's highest similarity to any term in the other
> network, averaged over both directions so each network counts equally regardless of size. To
> express each pairwise similarity as the fraction of the maximum possible similarity achievable
> above chance, we built a null distribution from 1,000 permutations per pair, each drawing two
> random term sets of the same sizes from the full set of GO Biological Process terms annotated
> genome-wide (org.Hs.eg.db; n = 12,588 terms) and recomputing the Wang/BMA similarity. The
> chance-corrected similarity was calculated as (observed similarity − mean null similarity) /
> (1 − mean null similarity). Empirical one-sided P values (the proportion of null-distribution
> values ≥ the observed similarity) were corrected for multiple testing across all 10 pairwise
> comparisons using the Benjamini-Hochberg procedure.

**Caveat**: the "full genome-wide GO:BP annotation (n = 12,588 terms)" sentence above is this
reconstruction's own choice for the permutation background, not something the legend specifies —
don't treat it as settled until the OCD–Depression discrepancy with the collaborator's
implementation is resolved (see above).

### Results, v3 (5 groups, BH-corrected)

| row | col | shared terms | raw similarity | expected (chance) | score | BH-adjusted P |
|---|---|---|---|---|---|---|
| dogCD | OCD | 0 | 0.185 | 0.207 | −0.03 | 0.838 |
| dogCD | OCD/dogCD | 5 | 0.432 | 0.225 | 0.27* | 0.00167 |
| dogCD | Depression | 5 | 0.396 | 0.221 | 0.22* | 0.00167 |
| dogCD | Schizophrenia | 5 | 0.352 | 0.222 | 0.17* | 0.00167 |
| OCD | OCD/dogCD | 0 | 0.120 | 0.224 | −0.13 | 1.000 |
| OCD | Depression | 7 | 0.221 | 0.220 | 0.00 | 0.679 |
| OCD | Schizophrenia | 0 | 0.118 | 0.222 | −0.13 | 1.000 |
| OCD/dogCD | Depression | 66 | 0.865 | 0.415 | 0.77* | 0.00167 |
| OCD/dogCD | Schizophrenia | 46 | 0.799 | 0.411 | 0.66* | 0.00167 |
| Depression | Schizophrenia | 50 | 0.724 | 0.442 | 0.51* | 0.00167 |

BH correction doesn't flip any significance calls here relative to raw P (the strong pairs stay
significant, OCD's pairs stay non-significant) — but it's now methodologically faithful to what
the legend specifies rather than an approximation.

Outputs: `data/results/05_network/pipeline_output/fig2/fig2e_similarity_heatmap_v3_full_legend.{R,png,pdf}`
and `fig2e_similarity_results_v3_full_legend.csv` — separate files again, nothing above was
overwritten.
