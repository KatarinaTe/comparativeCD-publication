# 03. QC

## Genotype concordance

We performed genotype concordance analyses across three different sets:

| Analysis | Number of dogs | QC analyses | 
|----------|-------------|-------------------|
| `1_High_LowPass` | 11 dogs | genotype concordance | 
| `2_Axiom_LowPass` | 5 dogs | genotype concordance | 
| `3_Downsampling` | 10 dogs | genotype concordance & depth | 


## 1_High vs LowPass:

11 dogs 
``` 
high_vs_lowpass.sh
``` 
Output files: genotype concordance for all pairs (we are reporting NON_REF_GENOTYPE_CONCORDANCE from these):
* Beskow_dup1.v2_concordance.vcf.genotype_concordance_summary_metrics
* Beskow.v2_concordance.vcf.genotype_concordance_summary_metrics
* Clarence.v2_concordance.vcf.genotype_concordance_summary_metrics
* Esme.v2_concordance.vcf.genotype_concordance_summary_metrics
* Finch.v2_concordance.vcf.genotype_concordance_summary_metrics
* Gus.v2_concordance.vcf.genotype_concordance_summary_metrics
* Jamie.v2_concordance.vcf.genotype_concordance_summary_metrics
* Kaylee_dup1.v2_concordance.vcf.genotype_concordance_summary_metrics
* Kaylee_dup2.v2_concordance.vcf.genotype_concordance_summary_metrics
* Kaylee_dup3.v2_concordance.vcf.genotype_concordance_summary_metrics
* Kaylee.v2_concordance.vcf.genotype_concordance_summary_metrics*
* Lily.v2_concordance.vcf.genotype_concordance_summary_metrics
* Lucky.v2_concordance.vcf.genotype_concordance_summary_metrics
* Rudy.v2_concordance.vcf.genotype_concordance_summary_metrics
* Tui.v2_concordance.vcf.genotype_concordance_summary_metrics

-----------------

## 2_Axiom vs LowPass:

5 dogs

``` 
axiom_vs_lowpass.sh

``` 
Output files: genotype concordance for all pairs (we are reporting NON_REF_GENOTYPE_CONCORDANCE from these):
* wholedog_gc_dogID10567_concordance.vcf.genotype_concordance_contingency_metrics 
* wholedog_gc_dogID10915_concordance.vcf.genotype_concordance_contingency_metrics 
* wholedog_gc_dogID12129_concordance.vcf.genotype_concordance_contingency_metrics 
* wholedog_gc_dogID2130_concordance.vcf.genotype_concordance_contingency_metrics 
* wholedog_gc_dogID7603concordance.vcf.genotype_concordance_contingency_metrics 
* wholedog_gc_dogNO_PAIR_concordance.vcf.genotype_concordance_contingency_metrics 
------------

## 3_Downsampling:

10 dogs

``` 
1_create_downsampled_bams.sh

2_mapping&imputing&depth_bams.sh

3_genotype_concordance.sh
``` 

Output files:
* **combined_knownsites.depth.txt** - depth for all bams
* **combined_concordance.txt** - genotype concordance summary for all samples

