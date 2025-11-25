process GBCMS_RUN {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/gbcms", mode: params.publish_dir_mode

    container "ghcr.io/msk-access/py-gbcms:2.1.2"

    input:
    tuple val(meta), path(bam), path(bai)
    path variants
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path("*.{vcf,maf}"), emit: counts
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = params.format ?: 'vcf'
    
    // Use per-sample suffix from meta, fallback to global params.suffix
    def suffix = meta.suffix ?: params.suffix
    
    // Construct filter arguments
    def filters = ""
    if (params.filter_duplicates)    filters += " --filter-duplicates"
    if (params.filter_secondary)     filters += " --filter-secondary"
    if (params.filter_supplementary) filters += " --filter-supplementary"
    if (params.filter_qc_failed)     filters += " --filter-qc-failed"
    if (params.filter_improper_pair) filters += " --filter-improper-pair"
    if (params.filter_indel)         filters += " --filter-indel"

    """
    gbcms run \\
        --variants ${variants} \\
        --bam ${prefix}:${bam} \\
        --fasta ${fasta} \\
        --output-dir . \\
        --format ${format} \\
        --suffix ${suffix} \\
        --threads ${task.cpus} \\
        --min-mapq ${params.min_mapq} \\
        --min-baseq ${params.min_baseq} \\
        ${filters} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gbcms: \$(python -c "import gbcms; print(gbcms.__version__)" 2>/dev/null || echo "2.0.0")
    END_VERSIONS
    """
}
