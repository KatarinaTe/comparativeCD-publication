#!/usr/bin/env python3
"""build_bam_dogid_keep_list.py — maps 01_mapping's BAM filenames (raw sequencing sampleIDs) to
the dogID convention used by 3_Merging_filtering's QC5/QC6 genotypes, producing a plink `--keep`
file for exactly the BAM-having dogs with a known genotype-ID mapping.

01_mapping's BAMs are named "<sampleID>.sorted.merged.MarkDups.BQSR.bam" (see
DOG10K_MERGE_MARKDUPS_BQSR_GVCF's output in 01_mapping/process_definitions.nf) — sampleID here is
the same raw sequencing sample identifier 3_Merging_filtering's RELABEL_AND_EXCLUDE_ROUND3 process
later relabels to "dogID". That relabel uses a manually-prepared crosswalk fam
(modi_DA_MERGED_GENCOVE_AXIOM_QC4C_forDogIDchange.fam) applied by row-position substitution, which
only carries the NEW dogIDs — no explicit sampleID<->dogID lookup table is published anywhere in
3_Merging_filtering's own output.

DA_MERGED_GENCOVE_AXIOM_QC3modi.fam (a raw, already-deposited data file — QC3 predates the
relabel) already carries both columns directly: `sampleID` and `IID` (the dogID), consistent
row-for-row with the relabel fam. So it doubles as the crosswalk this needs, with zero changes to
3_Merging_filtering's own code or published outputs required.
"""

import argparse
import re

parser = argparse.ArgumentParser(description="Map BAM-having sampleIDs to dogIDs via the QC3modi crosswalk fam.")
parser.add_argument("bam_dir_listing", help="A text file listing BAM filenames, one per line (e.g. `ls lowpass/bam/*.bam`)")
parser.add_argument("crosswalk_fam", help="DA_MERGED_GENCOVE_AXIOM_QC3modi.fam (has header: FID sampleID F M sex phe IID)")
parser.add_argument("out_keep_file", help="Output plink --keep file (FID IID, dogID-labeled)")
parser.add_argument("out_unmatched_file", help="Output list of BAM sampleIDs with no crosswalk match")
args = parser.parse_args()

BAM_SUFFIX = ".sorted.merged.MarkDups.BQSR.bam"

bam_sample_ids = set()
with open(args.bam_dir_listing) as f:
    for line in f:
        fname = line.strip().split("/")[-1]
        if not fname:
            continue
        if fname.endswith(BAM_SUFFIX):
            bam_sample_ids.add(fname[: -len(BAM_SUFFIX)])
        else:
            bam_sample_ids.add(re.sub(r"\.bam$", "", fname))

sample_to_dogid = {}
with open(args.crosswalk_fam) as f:
    header = f.readline().strip().split()
    sample_col = header.index("sampleID")
    iid_col = header.index("IID")
    for line in f:
        fields = line.strip().split()
        sample_to_dogid[fields[sample_col]] = fields[iid_col]

matched = []
unmatched = []
for sample_id in sorted(bam_sample_ids):
    dogid = sample_to_dogid.get(sample_id)
    if dogid is not None:
        matched.append(dogid)
    else:
        unmatched.append(sample_id)

with open(args.out_keep_file, "w") as f:
    for dogid in matched:
        f.write(f"{dogid} {dogid}\n")

with open(args.out_unmatched_file, "w") as f:
    for sample_id in unmatched:
        f.write(f"{sample_id}\n")

print(f"BAM-having sampleIDs: {len(bam_sample_ids)}")
print(f"Matched to a dogID via the crosswalk: {len(matched)}")
print(f"Unmatched (no crosswalk entry): {len(unmatched)}")
