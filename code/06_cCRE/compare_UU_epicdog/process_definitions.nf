// ============================================================================
// process_definitions.nf — 06. cCRE — UU vs. EpicDog element-count comparison
//
// Source: 06_cCRE/compare_UU_epicdog_elements.R (collaborator script, left untouched at its
// original location). bin/compare_UU_epicdog_elements.R is the wired copy -- inputs are
// staged under their real filenames (cerebellum_13_dense.bed, cerebrum_13_dense.bed,
// all_tissues_filtered.bed), so no CLI args were needed, only two real bugs fixed (see bin/
// header): a UCSC track-header line valr::read_bed() can't parse, and valr::read_bed()'s
// n_fields defaulting to 3, which left column 4 named X4 instead of name.
//
// New container: r-readr/r-dplyr/r-tidyr/r-ggplot2/r-valr, all conda-forge.
//
// BUILD_ALL_TISSUES_FILTERED (added 2026-09-23): live-builds the UU-side input
// (all_tissues_filtered.bed) from the 8 per-region UU BEDs instead of pointing at a static
// deposited copy -- source: 06_cCRE/build_all_tissues_filtered.sh (original, untouched),
// bin/build_all_tissues_filtered.sh is the wired copy. Reuses the gawk container already
// frozen for 01_mapping/2_cohort_genotyping -- no new build, this is plain awk/bash.
// ============================================================================

process BUILD_ALL_TISSUES_FILTERED {
    label 'process_single'
    container 'community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156'

    input:
    path(acg_bed,              arity: '1')
    path(cerebellum_bed,       arity: '1')
    path(frontal_lobe_bed,     arity: '1')
    path(hypothalamus_bed,     arity: '1')
    path(occipital_cortex_bed, arity: '1')
    path(striatum_bed,         arity: '1')
    path(temporal_cortex_bed,  arity: '1')
    path(thalamus_bed,         arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_all_tissues_filtered.sh
    """

    stub:
    """
    touch all_tissues_filtered.bed
    """

    output:
    path("all_tissues_filtered.bed"), emit: bed
}

process COMPARE_UU_EPICDOG_ELEMENTS {
    label 'process_single'
    container 'community.wave.seqera.io/library/compare_uu_epicdog:ebff580a4fa19c19'

    // No filename overrides needed -- Nextflow stages each input under its own original
    // basename by default, and the deposited/fetched files already ARE named
    // cerebellum_13_dense.bed / cerebrum_13_dense.bed / all_tissues_filtered.bed, which is
    // exactly what the untouched collaborator script's hardcoded read_dense_bed() calls
    // expect. See compare_UU_epicdog.nf for where those real filenames come from.
    input:
    path(epic_cerebellum_bed, arity: '1')
    path(epic_cerebrum_bed,   arity: '1')
    path(uu_all_tissues_bed,  arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    compare_UU_epicdog_elements.R
    """

    stub:
    """
    touch EPIC_comparison.txt
    touch Fig4d_corrected.pdf
    """

    output:
    path("EPIC_comparison.txt"),   emit: table
    path("Fig4d_corrected.pdf"),   emit: pdf
}
