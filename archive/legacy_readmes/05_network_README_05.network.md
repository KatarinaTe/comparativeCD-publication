# 05. network

We run network analyses using two github python environments:

https://github.com/ucsd-ccbb/NetColoc

https://github.com/sarah-n-wright/CrossSpeciesBMI

Scripts provided below are from these two githubs and edited to match our human/dog datasets.

Steps using iterations will create slightly different outputs between runs, output files are provided for these steps: i.e. colocalization enrichment incl control sets and HiDef systems map.

## A note on variable naming

CrossSpeciesBMI was originally written for a human-rat BMI comparison. Many variable, column, and
dictionary-key names throughout its scripts and the notebooks below still reflect that — most
visibly the `r`/`rat` naming (`NPS_r`, `rBMI`, `seed_bin_rat_BMI`, `dCCD` sometimes appearing as
`r...`) — and were never renamed when this project repurposed the code for a human-dog comparison.
Wherever you see `r`/`rat` in this code, read it as **dog**; `h`/`human` keeps its literal meaning.
The code confirms this itself, e.g. `2.1_OCD_dogCD_Network_Colocalization_260521.ipynb` cell 60:
`# all dog are called rat`.

## To run NetColoc/CrossSpeciesBMI, we need:

* **input seed genes** - for each trait, see below.
* **PCNet2.0** -  interactome_uuid='d73d6357-e87b-11ee-9621-005056ae23aa' (for PCNet2.0 Number of nodes: 19267 / Number of edges: 3852119)

### Precalculated matrices needed for propagation:
* **w_double_prime** - indiv_heats_OCDv4_CCDsept25.npy - can be rerun (see script) or use deposited 2.8Gb
* **w_prime** - w_prime_OCDv4_CCDsept25.npy - can be rerun (see script) or use deposited 2.8Gb

### Additonal tools:
* **NDEx** - https://www.ndexbio.org/
* **Cytoscape** - we used Cytoscape_v3.10.2

-------

## Network analyses run:

### Step 1. NetColoc:

The score cutoffs were set to  NPSd or NPSh >1.5 and NPShd > 3 for all three cross-species analyses.
(OBS! z-score = NPS)

#### 1a. Cross-species networks

| Analysis | Script name | Input: seed genes | Outputs: combined Z-scores | Colocalization enrichment | Systems map  |
|----------|-------------|-------------------|----------------------------|------------------------|-----------------|
| `OCD/dogCD` | 1.1_OCD_dogCD_NetColoc_analysis_260521.ipynb | OCDgenes_v4.txt / CCDgenes_Oct25.txt | ocd_ccd_zcomb_z12_251022.txt | netcoloc_enrichment_df_CCD_OCD_251022.csv  |  OCDv4_CCD_Oct25_systems_map.txt |
| `DEP/dogCD` | 1.1_OCD_dogCD_NetColoc_analysis_260521.ipynb (OBS! Change input) | human_DEP_noMHC_seed_genes_header.txt / CCDgenes_Oct25.txt | DEP_CCD_Oct25_zcomb_z12.txt | netcoloc_enrichment_df_CCD_DEP.csv  | DEP_CCD_Oct25_systems_map.txt | 
| `SCH/dogCD` | 1.1_OCD_dogCD_NetColoc_analysis_260521.ipynb (OBS! Change input) | human_schiz_seed_genes_header.txt / CCDgenes_Oct25.txt | SCH_CCD_zcomb_z12.txt |netcoloc_enrichment_df_CCD_SCH.csv  |  SCH_CCD_Oct25_systems_map.txt | 

#### 1b. Single species networks

NPS cutoff set to 3 for all single species analyses.

| Analysis | Script name | Input: seed genes | Outputs: Z-scores | Systems map | 
|----------|-------------|-------------------|------------------------------------|--------------|
| `OCD` | 1.2_single-species_OCD_NetColoc_analysis_260521.ipynb | OCDgenes_v4.txt | z_D1_OCD_z-scores_pcnet2_251022.csv | OCD_ONLY_systems_map.txt | 
| `DEP` | 1.2_single-species_OCD_NetColoc_analysis_260521.ipynb (OBS! Change input) | human_DEP_noMHC_seed_genes_header.txt | z_D1_DEP_z-scores_pcnet2.csv | DEP_ONLY_systems_map.txt |
| `SCH` | 1.2_single-species_OCD_NetColoc_analysis_260521.ipynb (OBS! Change input) | human_schiz_seed_genes_header.txt | z_D1_SCH_z-scores_pcnet2.csv | SCH_systems_map.txt |
| `dogCD` | 1.2_single-species_OCD_NetColoc_analysis_260521.ipynb (OBS! Change input) | CCDgenes_Oct25.txt| z_D2_CCD_Oct25_z-scores_pcnet2_251022.csv| CCD_ONLY_systems_map.txt |

---

### Step 2. CrossSpeciesBMI:
We need the seed genes and the output files from NetColoc, above.
Use new script in CrossSpeciesBMI, which produces the following outputs:


#### 2.1. Cross-species control analyses:
Run control analyses using seed genes, z-scores combined list, and this script:

| Analysis | Script name | Output: colocalization of control sets | Plots| 
|--------|-------------|-------------|-------------|
| `OCD/dogCD` | 2.1_OCD_dogCD_Network_Colocalization_260521.ipynb | controlanalyses_260521.csv |  Figure 1g-f | 

##### We need: control trait input seed genes:

| Control trait |  Input seed genes | Reference | 
|--------|-------------|-------------|
| `human height` | human_height2_seed_genes.txt |  | 
| `human RA` | human_RA_seed_genes.txt |  | 
| `human diabetes` | human_type2diabetes_seed_genes.txt |  | 
| `dog size` | DogHeight_seed_genes.txt |  | 
| `human height` | human_height2_seed_genes.txt |  | 
| `human depression` | human_DEP_noMHC_seed_genes.txt |  | 
| `human schizophrenia` | human_schiz_seed_genes.txt |  | 

---
#### 2.2. Cross-species networks

| Analysis | Script name | Outputs: GO-term enrichment | MPO-term enrichment  |   
|--------|-------------|-------------|-------------|
| `OCD/dogCD` | 2.2_OCD_dogCD_Systems_Map_260521.ipynb | Compulsive_hierachy_full_GO_enrichment_251023.tsv | 251023_KT_hierarchy_full_MGD_enrichment_results.tsv| 
| `DEP/dogCD` | 2.2_OCD_dogCD_Systems_Map_260521.ipynb (OBS! Change input) | DEP_CCD_hierachy_full_GO_enrichment.tsv | DEP_CCD_hierarchy_full_MGD_enrichment_results.tsv| 
| `SCH/dogCD` | 2.2_OCD_dogCD_Systems_Map_260521.ipynb (OBS! Change input) | SCH_CCD_hierachy_full_GO_enrichment.tsv | SCH_CCD_hierarchy_full_MGD_enrichment_results.tsv| 



#### 2.3. Single species networks

| Analysis | Script name | Outputs: GO-term enrichment | MPO-term enrichment  |   
|--------|-------------|-------------|-------------|
| `OCD` | 2.3_OCD_ONLY_Systems_Map_260525.ipynb | OCD_ONLY_hierachy_full_GO_enrichment.tsv | OCD_ONLY_hierarchy_full_MGD_enrichment_results.tsv| 
| `DEP` | 2.3_OCD_ONLY_Systems_Map_260525.ipynb (OBS! Change input) | DEP_ONLY_hierarchy_full_GO_enrichment.tsv | DEP_ONLY_hierarchy_full_MGD_enrichment_results.tsv| 
| `SCH` | 2.3_OCD_ONLY_Systems_Map_260525.ipynb (OBS! Change input) | SCH_ONLY_hierachy_full_GO_enrichment.tsv | SCH_ONLY_hierarchy_full_MGD_enrichment_results.tsv| 
| `dogCD` | 2.3_OCD_ONLY_Systems_Map_260525.ipynb (OBS! Change input) | CCD_ONLY_hierachy_full_GO_enrichment.tsv | CCD_ONLY_hierarchy_full_MGD_enrichment_results.tsv| 

---


#### 3. Cytoscape - network (systems map and gene communities) visualization

We used Cytoscape Version: 3.10.2 for visalization purposes, all system maps, HiDef and cosine similarities for subcommunities for all cross-species and single-species networks are all incluced in here: 

* **CCD_Oct25_system250522_publication.cys** 

The rest of Figure editing was done in Adobe illustrator v.26.0.2


#### Original files from ndexbio.org (used as inputs for Cytoscape) are called: (located in data/05_network/cytoscape):
* **prefix_NetColoc_subgraph_CosSim95.cx** - cosine similarity
* **prefix_NetColoc_subgraph.cx** - HiDef
* **prefix_systems_map.cx** - hierchial systems map

| Analysis | Prefix | 
|--------|-------------|
| `OCD/dogCD` | OCDv4_CCD_Oct25 | 
| `depression/dogCD` | DEP_CCD_Oct25 |
| `schizophrenia/dogCD` | SCH_CCD_Oct25 |  
| `OCD` | OCD | 
| `depression` | DEP |
| `schizophrenia` | SCH |  
| `dogCD` | CCD| 



---
## /Figures

| Figure | Script name | Source data | High resolution image |
|--------|-------------|-------------|-----------------------|
| Figure 1e | plot_netcoloc_zscores.R  | data/Figures |   | 
| Figure 1g-f | 1_OCD_dogCD_Network_Colocalization_260521.ipynb |  |   | 
