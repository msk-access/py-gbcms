#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.variants) { ch_variants = file(params.variants) } else { exit 1, 'Variants file not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Reference FASTA not specified!' }

/*
========================================================================================
    IMPORT LOCAL MODULES/WORKFLOWS
========================================================================================
*/

include { GBCMS } from './workflows/gbcms'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    Channel
        .fromPath(params.input)
        .splitCsv(header:true, sep:',')
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            
            // Per-sample suffix: use row.suffix if present, otherwise use global params.suffix
            meta.suffix = row.containsKey('suffix') && row.suffix ? row.suffix : params.suffix
            
            def bam = file(row.bam, checkIfExists: true)
            
            // Handle BAI: if provided use it, otherwise auto-discover with both naming conventions
            def bai
            if (row.bai) {
                bai = file(row.bai, checkIfExists: true)
            } else {
                // Try both common BAI naming conventions: .bam.bai and .bai
                def bai_path1 = "${row.bam}.bai"
                def bai_path2 = row.bam.replaceAll(/\.bam$/, '.bai')
                def bai1 = file(bai_path1)
                def bai2 = file(bai_path2)
                
                if (bai1.exists()) {
                    bai = bai1
                } else if (bai2.exists()) {
                    bai = bai2
                } else {
                    error "BAI index not found for ${row.bam}. Searched: ${bai_path1}, ${bai_path2}"
                }
            }
            
            return [ meta, bam, bai ]
        }
        .set { ch_samplesheet }
    
    // Prepare other inputs
    ch_variants_file = file(params.variants)
    
    // FASTA and FAI
    ch_fasta_file = file(params.fasta)
    ch_fai_file = file("${params.fasta}.fai")
    ch_fasta_tuple = [ ch_fasta_file, ch_fai_file ]

    //
    // WORKFLOW: Run main workflow
    //
    GBCMS (
        ch_samplesheet,
        ch_variants_file,
        ch_fasta_tuple
    )
}
