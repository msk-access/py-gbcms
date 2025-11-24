/*
========================================================================================
    IMPORT MODULES / SUBWORKFLOWS
========================================================================================
*/

include { GBCMS_RUN } from '../modules/local/gbcms/run/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow GBCMS {
    take:
    ch_samplesheet // channel: [ val(meta), [ bam, bai ] ]
    ch_variants    // channel: [ path(variants) ]
    ch_fasta       // channel: [ path(fasta), path(fai) ]

    main:
    
    ch_versions = Channel.empty()

    //
    // MODULE: Run gbcms
    //
    GBCMS_RUN (
        ch_samplesheet,
        ch_variants,
        ch_fasta
    )
    ch_versions = ch_versions.mix(GBCMS_RUN.out.versions)

    emit:
    counts   = GBCMS_RUN.out.counts
    versions = ch_versions
}
