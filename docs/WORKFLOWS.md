# Usage Patterns

`py-gbcms` can be used in two ways depending on your needs:

## 🔧 Standalone CLI (Single/Few Samples)

**Best for:**
- Processing 1-10 samples
- Quick ad-hoc analysis
- Local development and testing
- Direct control over parameters

**Pros:**
- ✅ Simple - just one command
- ✅ Fast to set up
- ✅ Full control over threading
- ✅ Works anywhere (local, server, container)

**Cons:**
- ❌ Manual parallelization for multiple samples
- ❌ No automatic resource management
- ❌ Requires manual error handling

**Example:**
```bash
gbcms run \
    --variants variants.vcf \
    --bam sample1.bam \
    --fasta reference.fa \
    --output-dir results/
```

**Learn more:** [CLI Quick Start](quick-start.md)

---

## 🔄 Nextflow Workflow (Many Samples, HPC)

**Best for:**
- Processing 10+ samples
- HPC/SLURM cluster environments
- Reproducible pipelines
- Automated retry and error handling

**Pros:**
- ✅ Automatic parallelization across samples
- ✅ Smart resource management
- ✅ Built-in retry logic
- ✅ Resume failed runs
- ✅ Portable (Docker/Singularity)

**Cons:**
- ❌ Requires Nextflow installation
- ❌ More setup (samplesheet, config)
- ❌ Learning curve for Nextflow DSL

**Example:**
```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    -profile slurm
```

**Learn more:** [Nextflow Workflow Guide](NEXTFLOW.md)

---

## 📊 Quick Comparison

| Feature | CLI | Nextflow |
|---------|-----|----------|
| Setup complexity | Low | Medium |
| Best for # samples | 1-10 | 10+ |
| Parallelization | Manual | Automatic |
| Resource management | Manual | Automatic |
| HPC integration | Manual | Built-in |
| Resume failed jobs | No | Yes |
| Reproducibility | Good | Excellent |

## Which Should I Use?

### Use **CLI** if you:
- Have a few samples to process
- Want quick results
- Are working locally or on a single server
- Need full manual control

### Use **Nextflow** if you:
- Have many samples (10+)
- Are on an HPC cluster (SLURM, PBS, etc.)
- Need reproducible pipelines
- Want automatic parallelization and error handling
- Plan to re-run the analysis multiple times
