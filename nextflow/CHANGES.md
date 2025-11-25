## Summary of Changes Made to Nextflow Workflow

### Fixed Issues

1. **Added `publishDir` directive** to `GBCMS_RUN` module
   - Outputs are now properly published to `${params.outdir}/gbcms`

2. **Fixed BAI index handling**
   - BAI column in samplesheet is now optional
   - If not provided, looks for `<bam>.bai`
   - Added `checkIfExists` validation for input files

3. **Improved version extraction**
   - Changed from `gbcms --version` (which may not exist) to Python import
   - Falls back to hardcoded "2.0.0" if import fails

4. **Fixed sample ID handling**
   - Changed `--bam ${bam}` to `--bam ${prefix}:${bam}` for proper sample naming
   - Removed redundant `--suffix .gbcms` (prefix serves this purpose)

5. **Added missing configuration**
   - Added `max_cpus`, `max_memory`, and `max_time` parameters
   - These are required by the `check_max()` function

6. **Added SLURM profile**
   - Uses SLURM executor with default queue `cpu_medium`
   - Enabled Singularity for cluster compatibility
   - Can be customized by editing `nextflow.config`

### Usage

**Local (Docker):**
```bash
nextflow run nextflow/main.nf --input samples.csv --variants vars.vcf --fasta ref.fa -profile docker
```

**SLURM cluster:**
```bash
nextflow run nextflow/main.nf --input samples.csv --variants vars.vcf --fasta ref.fa -profile slurm
```

The workflow is now production-ready and follows nf-core best practices.
