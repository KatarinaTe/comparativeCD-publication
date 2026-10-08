import pandas as pd
import argparse

# Hardcoded bim file path
bim_path = "/path/to/CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.bim"

# Set up argument parser
parser = argparse.ArgumentParser(description="Harmonize sumstats alleles and BETA to reference BIM file.")
parser.add_argument("sumstats_in", help="Input sumstats file")
parser.add_argument("sumstats_out", help="Output harmonized sumstats file")
args = parser.parse_args()

# Load sumstats
sumstats = pd.read_csv(args.sumstats_in, sep="\s+", dtype=str)
# Load bim
bim = pd.read_csv(bim_path, sep="\s+", header=None, names=["CHR", "SNP", "CM", "POS", "A1_ref", "A2_ref"], dtype=str)

# Merge on SNP ID
merged = sumstats.merge(bim[["SNP", "A1_ref", "A2_ref"]], on="SNP", how="inner")

# Track swap status
def harmonize(row):
    if row["A1"] == row["A1_ref"] and row["A2"] == row["A2_ref"]:
        return row["A1"], row["A2"], row["BETA"], False
    elif row["A1"] == row["A2_ref"] and row["A2"] == row["A1_ref"]:
        return row["A2"], row["A1"], str(-float(row["BETA"])), True
    else:
        return pd.NA, pd.NA, pd.NA, pd.NA

merged[["A1_harmonized", "A2_harmonized", "BETA_harmonized", "swapped"]] = merged.apply(
    harmonize, axis=1, result_type="expand"
)

total_snps = merged.shape[0]
swapped_count = merged["swapped"].sum(skipna=True)
dropped_count = merged["A1_harmonized"].isna().sum()

harmonized = merged.dropna(subset=["A1_harmonized"])
harmonized["A1"] = harmonized["A1_harmonized"]
harmonized["A2"] = harmonized["A2_harmonized"]
harmonized["BETA"] = harmonized["BETA_harmonized"]
harmonized = harmonized[sumstats.columns]

harmonized.to_csv(args.sumstats_out, sep="\t", index=False)

# Save list of swapped SNPs
swapped_snps = merged.loc[merged["swapped"] == True, "SNP"]
swapped_snps.to_csv("swapped_snps.txt", index=False, header=False)

print(f"Total SNPs processed: {total_snps}")
print(f"Alleles swapped (and BETA flipped): {int(swapped_count)}")
print(f"SNPs dropped due to mismatched alleles: {int(dropped_count)}")
print(f"SNPs retained in output: {harmonized.shape[0]}")
print(f"List of swapped SNPs saved to swapped_snps.txt")
