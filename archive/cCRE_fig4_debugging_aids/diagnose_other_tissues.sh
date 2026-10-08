#!/bin/bash
set -euo pipefail

DIR="data/06_cCRE/epicdog_chromatin_states_bed_canfam4"

echo "--- per-tissue active-state bp and chrom count ---"
for t in colon kidney liver lung mammary_gland ovary pancreas spleen stomach; do
  f="${DIR}/${t}_13_dense.bed"
  n_chr=$(tail -n +2 "$f" | cut -f1 | sort -u | wc -l)
  bp=$(tail -n +2 "$f" | awk 'BEGIN{FS="\t"} $4 ~ /^(7|8|9|11|12|13)$/ {sum+=$3-$2} END{print sum+0}')
  echo -e "${t}\tunique_chroms=${n_chr}\tactive_bp=${bp}"
done

echo ""
echo "--- chrom name samples per file (first 3 unique) ---"
for t in cerebellum cerebrum colon kidney liver lung mammary_gland ovary pancreas spleen stomach; do
  f="${DIR}/${t}_13_dense.bed"
  echo -n "${t}: "
  tail -n +2 "$f" | cut -f1 | sort -u | head -3 | tr "\n" " "
  echo
done
