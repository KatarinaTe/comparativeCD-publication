#!/usr/bin/env python3
"""build_st5_st6_region_ids.py -- standalone, non-pipeline script (not wired into any
Nextflow process). Verifies and reproduces the region/gene counts in the manuscript's
Results text ("27 regions... significant, and 28 additional regions... suggestive"; "46"/
"41" genes; the chromosome-1 region with 10 genes) directly from Supplementary Table 5, and
builds the columns that were pasted into ST5 and ST6 to make those counts traceable.

Region-merging rule (confirmed against the manuscript, 2026-09-29): two of ST5's rows are
the same genomic region if their "Wide region (+/- 100kb)" intervals overlap on the same
chromosome. This is the only merge rule that reproduces the real, correct "27 significant
regions" figure exactly (a naive per-row count over-counts, since the same physical locus
can be hit independently by multiple survey items/factors). A merged region is classified
significant if ANY of its constituent rows is significant, suggestive only if none are --
so a suggestive-only association (e.g. CNGB1) that happens to share a merged region with a
different gene's significant hit still traces to that region's sig_region_number, even
though CNGB1's own classification (used for the 46/41 gene tally) stays "suggestive".

Gene assignment is unaffected by this region-level merging: genes come from each row's own
credible-set interval (ST5 column S, "Genes (NetColoc)"), never from the wider merged
region -- this is why the chr1 region (merged_region_id 3) has exactly 10 genes (its own
CS), not the much larger gene set that would fall within its full +/-100kb clump window.

Outputs (not deposited/committed by this script -- see its own run instructions):
  ST5_regions_with_merged_ids.csv/.xlsx   -- ST5 + merged_region_id/range/sig_region_number/
                                              sugg_region_number columns
  ST5_gene_to_region_trace.csv/.xlsx      -- one row per gene (of the 87 in ST5 column S),
                                              its classification, region trace, and whether
                                              it's in the deposited 83-gene dogCD seed list
  ST6_with_region_numbers.csv/.xlsx       -- ST6 (Human and dog genes taken forward for
                                              network analysis) + dogCD sig_region_number/
                                              dogCD sugg_region_number columns, 'NA' filled
                                              for non-applicable cells (e.g. human genes)

Usage:
  python3 build_st5_st6_region_ids.py <path to NATURE_SupplementaryTables.xlsx> \\
      <path to CCDgenes_Oct25.txt> <output directory>
"""
import sys
import openpyxl
import pandas as pd

PLACEHOLDER = {"no genes", "none", "na", "n/a", ""}


def load_seed_genes(path):
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    return set(lines[1:])  # skip header row ("gene")


def read_sheet(xlsx_path, sheet_name):
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb[sheet_name]
    rows = list(ws.iter_rows(values_only=True))
    return rows


def build_st5(xlsx_path, seed_genes, out_dir):
    rows = read_sheet(xlsx_path, "ST5")
    header = list(rows[0])
    newheader = [f"{h}...{i+1}" if header.count(h) > 1 else h for i, h in enumerate(header)]
    df = pd.DataFrame(rows[1:], columns=newheader)
    col_genes = newheader[18]  # column S, "Genes (NetColoc)"

    def parse_range(s):
        chrom, rest = s.split(":")
        start, end = rest.split("..")
        return chrom, float(start), float(end)

    df[["chrom", "wstart", "wend"]] = df["Wide region (+/- 100kb)"].apply(lambda s: pd.Series(parse_range(s)))
    df = df.sort_values(["chrom", "wstart"]).reset_index(drop=True)

    # ---- merge overlapping wide-regions within each chromosome ----
    merged = []
    cur = None
    for i, r in df.iterrows():
        if cur is None or r["chrom"] != cur["chrom"] or r["wstart"] > cur["wend"]:
            if cur is not None:
                merged.append(cur)
            cur = {"chrom": r["chrom"], "wstart": r["wstart"], "wend": r["wend"], "row_idx": [i]}
        else:
            cur["wend"] = max(cur["wend"], r["wend"])
            cur["row_idx"].append(i)
    if cur is not None:
        merged.append(cur)

    df["merged_region_id"] = 0
    df["merged_region_range"] = ""
    df["sig_region_number"] = pd.NA
    df["sugg_region_number"] = pd.NA
    sig_counter = sugg_counter = 0
    for m_i, m in enumerate(merged, start=1):
        region_range = f"{m['chrom']}:{int(m['wstart'])}-{int(m['wend'])}"
        is_sig = any(df.loc[i, "threshold"] == "significant" for i in m["row_idx"])
        if is_sig:
            sig_counter += 1
            sig_num, sugg_num = sig_counter, pd.NA
        else:
            sugg_counter += 1
            sig_num, sugg_num = pd.NA, sugg_counter
        for i in m["row_idx"]:
            df.loc[i, "merged_region_id"] = m_i
            df.loc[i, "merged_region_range"] = region_range
            df.loc[i, "sig_region_number"] = sig_num
            df.loc[i, "sugg_region_number"] = sugg_num

    print(f"Merged regions: {len(merged)} (significant 1-{sig_counter}, suggestive-only 1-{sugg_counter})")

    df["Region number"] = df["Region number"].astype(float)
    df = df.sort_values("Region number").drop(columns=["chrom", "wstart", "wend"]).reset_index(drop=True)

    out_cols = list(newheader) + ["merged_region_id", "merged_region_range", "sig_region_number", "sugg_region_number"]
    df_out = df[out_cols]
    df_out.to_csv(f"{out_dir}/ST5_regions_with_merged_ids.csv", index=False)
    with pd.ExcelWriter(f"{out_dir}/ST5_regions_with_merged_ids.xlsx", engine="openpyxl") as w:
        df_out.to_excel(w, sheet_name="ST5_with_region_ids", index=False)

    # ---- gene-to-region trace: classification uses the gene's OWN row threshold (this is
    # ---- what reproduces the manuscript's exact 46 significant / 41 suggestive gene counts),
    # ---- but region numbers point to whatever the gene's row actually carries -- these can
    # ---- differ (see CNGB1: classification stays "suggestive", but its merged region also
    # ---- contains a different gene's significant hit, so it has a sig_region_number too).
    gene_rows = {}
    for _, r in df.iterrows():
        v = r[col_genes]
        if pd.isna(v):
            continue
        for g in str(v).split(","):
            g = g.strip()
            if not g or g.lower() in PLACEHOLDER:
                continue
            e = gene_rows.setdefault(g, {"regions": set(), "own_sig": False, "own_sugg": False,
                                          "sig_nums": set(), "sugg_nums": set()})
            e["regions"].add(r["merged_region_id"])
            if r["threshold"] == "significant":
                e["own_sig"] = True
            else:
                e["own_sugg"] = True
            if pd.notna(r["sig_region_number"]):
                e["sig_nums"].add(int(r["sig_region_number"]))
            if pd.notna(r["sugg_region_number"]):
                e["sugg_nums"].add(int(r["sugg_region_number"]))

    gene_table = []
    for g, e in sorted(gene_rows.items()):
        classification = "significant" if e["own_sig"] else "suggestive"
        shared_locus = e["own_sugg"] and not e["own_sig"] and len(e["sig_nums"]) > 0
        gene_table.append({
            "gene": g,
            "classification": classification,
            "merged_region_ids": ";".join(str(x) for x in sorted(e["regions"])),
            "sig_region_number": ";".join(str(n) for n in sorted(e["sig_nums"])),
            "sugg_region_number": ";".join(str(n) for n in sorted(e["sugg_nums"])),
            "note": "suggestive hit shares a merged region with a different gene's significant hit" if shared_locus else "",
            "in_dogCD_seed_list_83": g in seed_genes,
        })
    gene_df = pd.DataFrame(gene_table)
    print(f"Gene-to-region trace: {len(gene_df)} genes "
          f"({(gene_df['classification']=='significant').sum()} significant, "
          f"{(gene_df['classification']=='suggestive').sum()} suggestive)")
    gene_df.to_csv(f"{out_dir}/ST5_gene_to_region_trace.csv", index=False)
    with pd.ExcelWriter(f"{out_dir}/ST5_gene_to_region_trace.xlsx", engine="openpyxl") as w:
        gene_df.to_excel(w, sheet_name="gene_to_region_trace", index=False)

    return gene_df


def build_st6(xlsx_path, gene_df, out_dir):
    rows = read_sheet(xlsx_path, "ST6")
    header = list(rows[1])
    data = [r for r in rows[2:] if r[0] is not None]
    df = pd.DataFrame(data, columns=header)
    df = df.drop(columns=[c for c in df.columns if c is None])
    df = df[df["Species"].isin(["Dog", "Human"])].reset_index(drop=True)  # drop footnote rows

    lookup = gene_df[["gene", "sig_region_number", "sugg_region_number"]].rename(columns={
        "gene": "Gene*",
        "sig_region_number": "dogCD sig_region_number",
        "sugg_region_number": "dogCD sugg_region_number",
    })
    merged = df.merge(lookup, on="Gene*", how="left")
    for col in ["dogCD GWAS region", "dogCD sig_region_number", "dogCD sugg_region_number"]:
        merged[col] = merged[col].apply(
            lambda v: "NA" if (pd.isna(v) or v == "")
            else (str(int(v)) if isinstance(v, (int, float)) and float(v).is_integer() else str(v))
        )

    merged.to_csv(f"{out_dir}/ST6_with_region_numbers.csv", index=False)
    with pd.ExcelWriter(f"{out_dir}/ST6_with_region_numbers.xlsx", engine="openpyxl") as w:
        merged.to_excel(w, sheet_name="ST6_with_regions", index=False)
    print(f"ST6 with region numbers: {merged.shape[0]} rows")


if __name__ == "__main__":
    xlsx_path, seed_path, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    seed_genes = load_seed_genes(seed_path)
    gene_df = build_st5(xlsx_path, seed_genes, out_dir)
    build_st6(xlsx_path, gene_df, out_dir)
