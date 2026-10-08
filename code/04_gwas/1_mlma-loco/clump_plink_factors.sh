plink --dog --bfile CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6  --out CCDF1_clump250 --clump CCDF1_QC6_LOCO.loco.mlma --clump-field p --clump-p1 0.000001 --clump-p2 0.00001 --clump-kb 250 --clump-r2 0.50 --clump-range genes6_UU_Cfam_GSD_1.0_ROSY.txt 
plink --dog --bfile CCDF2_DA_MERGED_GENCOVE_AXIOM_QC6  --out CCDF2_clump250 --clump CCDF2_QC6_LOCO.loco.mlma --clump-field p --clump-p1 0.000001 --clump-p2 0.00001 --clump-kb 250 --clump-r2 0.50 --clump-range genes6_UU_Cfam_GSD_1.0_ROSY.txt 
plink --dog --bfile CCDF3_DA_MERGED_GENCOVE_AXIOM_QC6  --out CCDF3_clump250 --clump CCDF3_QC6_LOCO.loco.mlma --clump-field p --clump-p1 0.000001 --clump-p2 0.00001 --clump-kb 250 --clump-r2 0.50 --clump-range genes6_UU_Cfam_GSD_1.0_ROSY.txt 

