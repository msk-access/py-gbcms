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
            
            def bam = file(row.bam, checkIfExists: true)
            
            // Handle BAI: if provided use it, otherwise auto-discover and validate
            def bai
            if (row.bai) {
                bai = file(row.bai, checkIfExists: true)
            } else {
                def bai_path = "${row.bam}.bai"
                bai = file(bai_path)
                if (!bai.exists()) {
                    error "BAI index not found for ${row.bam}. Expected: ${bai_path}"
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
