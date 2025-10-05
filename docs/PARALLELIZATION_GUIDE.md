# Parallelization Guide: joblib vs Ray

Complete guide to choosing and using the right parallelization backend for your gbcms workload.

## 🎯 Quick Reference
| Backend | Best For | Scale | Setup | Use Case |
|---------|----------|-------|-------|----------|
| **joblib** | Local multi-core | 1-64 cores | Simple | Most workloads |
| **Ray** | Multi-node clusters | 64-1000+ cores | Complex | Massive datasets |

### When to Use Each Backend

#### **Override Options**

**Generic Counting** (`--generic-counting`):
- **Purpose**: Force all variants to use generic counting algorithm
- **When to use**: When you need maximum compatibility or debugging complex variants
- **Trade-off**: Slower than smart strategy but works for all variant types
- **Use case**: Troubleshooting or when specialized algorithms fail

```bash
# Force generic counting for all variants
gbcms count run \
    --fasta reference.fa \
    --bam samples.txt \
    --vcf variants.vcf \
    --generic-counting \
    --thread 16
```

### joblib (Recommended for most users)

**✅ Use joblib when:**
- Processing <1M variants
- Single machine with multiple cores (2-64 cores)
- Local workstation or server
- Quick turnaround needed
- Memory usage is manageable

**📈 Performance:**
- **8 cores**: ~100x speedup vs single-threaded
- **32 cores**: ~200x speedup vs single-threaded
- **Memory**: ~1.2-2x baseline

**🔧 Backend Options:**
- `loky` (default): Best for CPU-intensive tasks
- `threading`: Best for I/O-intensive tasks
- `multiprocessing`: True parallelism, higher memory

### Ray (For massive scale)

**✅ Use Ray when:**
- Processing >1M variants
- Multi-node cluster available
- Long-running jobs (hours/days)
- Need fault tolerance
- Memory distributed across nodes

**📈 Performance:**
- **32 workers**: ~270x speedup vs single-threaded
- **128 workers (4 nodes)**: ~675x speedup vs single-threaded
- **Memory**: Distributed across cluster

## 🚀 joblib Usage Guide

### Installation

joblib is included with core gbcms installation - no extra setup needed.

### CLI Usage

```bash
# Use joblib (default backend)
gbcms count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output counts.txt \
    --thread 16 \
    --backend joblib

# Use threading for I/O-heavy workloads
gbcms count run \
    --thread 8 \
    --backend threading \
    ...
```

### Python API Usage

#### 1. Simple Parallel Processing

```python
from gbcms.parallel import parallel_map

def count_variant(variant):
    # Your counting logic here
    return counts

# Process variants in parallel
results = parallel_map(
    func=count_variant,
    items=variants,
    n_jobs=8,  # Use 8 cores
    backend='loky',  # or 'threading', 'multiprocessing'
    description="Counting variants",
    show_progress=True,
)
```

#### 2. Batch Processing (Memory Efficient)

```python
from gbcms.parallel import BatchProcessor

# Process in batches for better memory usage
processor = BatchProcessor(
    batch_size=1000,  # 1000 variants per batch
    n_jobs=8,
    backend='loky',
)

def process_batch(variant_batch):
    # Process entire batch at once
    return [count_variant(v) for v in variant_batch]

results = processor.process_batches(
    func=process_batch,
    items=all_variants,
    description="Processing variant batches",
)

processor.shutdown()
```

#### 3. Multiple Arguments (Starmap)

```python
from gbcms.parallel import parallel_starmap

def count_variant_with_config(variant, config, reference):
    # Function that needs multiple arguments
    return counts

# Create argument tuples
args_list = [(v, config, reference) for v in variants]

results = parallel_starmap(
    func=count_variant_with_config,
    items=args_list,
    n_jobs=8,
    description="Counting with config",
)
```

## 🌐 Ray Usage Guide

### Installation

```bash
# Install gbcms with Ray support
uv pip install "gbcms[ray]"

# Or with all features
uv pip install "gbcms[all]"
```

### CLI Usage

```bash
# Use Ray for distributed processing
gbcms count run \
    --fasta reference.fa \
    --bam-fof bam_files.txt \
    --vcf variants.vcf \
    --output counts.txt \
    --thread 32 \
    --backend ray \
    --use-ray

# Connect to existing Ray cluster
RAY_ADDRESS='ray://cluster-head:10001' gbcms count run \
    --thread 64 \
    --backend ray \
    --use-ray \
    ...
```

### Python API Usage

#### 1. Basic Distributed Processing

```python
from gbcms.parallel import ParallelProcessor

# Initialize with Ray (can exceed local CPU count)
processor = ParallelProcessor(
    n_jobs=32,  # Can be more than local CPUs
    backend='ray',
    use_ray=True,
)

# Automatically distributes across available workers
results = processor.map(
    func=count_variant,
    items=variants,
    description="Distributed counting",
)

processor.shutdown()
```

#### 2. Ray Actors (Stateful Workers)

```python
from gbcms.parallel import create_ray_actors, distribute_work_to_actors

# Create Ray actors (one per worker)
actors = create_ray_actors(
    n_actors=16,
    config_dict=config.dict(),
)

# Distribute work to actors
results = distribute_work_to_actors(
    actors=actors,
    work_items=variant_blocks,
    description="Processing with actors",
)
```

#### 3. Multi-Node Cluster Setup

```python
import ray

# Option 1: Connect to existing cluster
ray.init(address='ray://cluster-head:10001')

# Option 2: Start local cluster
ray.init(num_cpus=64, num_gpus=4)

# Use as normal
processor = ParallelProcessor(use_ray=True)
results = processor.map(count_variant, variants)
```

### Ray Dashboard

Monitor your Ray cluster:

```bash
# Start Ray head node with dashboard
ray start --head --dashboard-host=0.0.0.0

# Access dashboard at http://localhost:8265
```

## 🎯 High-Sample-Count Scenarios (100s-1000s of Samples)

Processing hundreds to thousands of samples is one of the most common gbcms use cases. This section covers specific optimizations for these workloads.

### When You're Processing 100s-1000s of Samples

**✅ This scenario applies when:**
- **100-1000+ BAM files** per analysis
- **Multiple variant files** (VCF/MAF) being processed together
- **Clinical sequencing** with many patient samples
- **Population studies** with large cohorts
- **Batch processing** of multiple experiments

**📊 Typical Performance:**
- **100 samples**: 10-30 minutes with joblib (16-32 threads)
- **500 samples**: 30-90 minutes with joblib (32 threads) or Ray (64+ workers)
- **1000+ samples**: 1-3 hours with Ray distributed processing

### joblib Strategy for High Sample Counts

#### 1. **Optimal Thread Configuration**

```bash
# For 100-500 samples: Use most cores but leave some for OS
gbcms count run \
    --fasta reference.fa \
    --bam-fof sample_list.txt \
    --vcf variants.vcf \
    --output results.txt \
    --thread 24 \
    --backend joblib \
    --max-block-size 5000
```

#### 2. **Memory-Efficient Processing**

```python
from gbcms.parallel import BatchProcessor

# Process samples in smaller batches to manage memory
processor = BatchProcessor(
    batch_size=50,  # 50 samples per batch
    n_jobs=16,
    backend='loky',
)

def process_sample_batch(sample_batch):
    results = []
    for sample in sample_batch:
        # Load BAM and process variants for this sample
        bam_file = sample_bam_mapping[sample]
        alignments = load_alignments(bam_file, variants)
        counts = count_variants_batch(variants, alignments, sample)
        results.append(counts)
    return results

# Process all samples in batches
all_results = []
for i in range(0, len(all_samples), 50):
    batch = all_samples[i:i+50]
    batch_results = processor.process_batches(process_sample_batch, [batch])
    all_results.extend(batch_results)
```

#### 3. **BAM File-of-Files Optimization**

```bash
# Create efficient sample list
cat > high_throughput_samples.txt << EOF
sample001	sample001.bam
sample002	sample002.bam
# ... 100s-1000s of samples
sample999	sample999.bam
EOF

# Process with optimized settings
gbcms count run \
    --fasta hg38.fa \
    --bam-fof high_throughput_samples.txt \
    --vcf panel.vcf \
    --thread 32 \
    --max-block-size 10000 \
    --backend joblib
```

### Ray Strategy for High Sample Counts

#### 1. **Distributed Sample Processing**

```python
from gbcms.parallel import ParallelProcessor

# For 1000+ samples, use Ray for distributed processing
processor = ParallelProcessor(
    n_jobs=64,  # Can exceed local CPU count
    backend='ray',
    use_ray=True,
)

# Distribute samples across cluster
sample_batches = [samples[i:i+100] for i in range(0, len(samples), 100)]

results = processor.map(
    func=process_sample_batch,
    items=sample_batches,
    description="Processing 1000+ samples across cluster",
)
```

#### 2. **Multi-Node Sample Distribution**

```bash
# For 1000+ samples across multiple nodes
# Node 1: 32 cores, Node 2: 32 cores, Node 3: 32 cores = 96 total workers

# Start Ray cluster across nodes
# Node 1 (head):
ray start --head --num-cpus=32 --dashboard-host=0.0.0.0

# Node 2:
ray start --address=ray://node1:10001 --num-cpus=32

# Node 3:
ray start --address=ray://node1:10001 --num-cpus=32

# Process with all 96 workers
gbcms count run \
    --bam-fof all_samples.txt \
    --thread 96 \
    --backend ray \
    --use-ray
```

### Memory Management for High Sample Counts

#### 1. **Per-Sample Memory Estimation**

```python
# Estimate memory per sample
def estimate_sample_memory(sample_name, variant_count):
    # Rough estimate: 100MB per sample for 10K variants
    base_memory = 100 * 1024 * 1024  # 100MB base
    per_variant_memory = 10 * 1024    # 10KB per variant
    return base_memory + (variant_count * per_variant_memory)

# For 500 samples × 10K variants = ~5GB total
total_memory = sum(estimate_sample_memory(s, 10000) for s in samples)
```

#### 2. **Memory-Optimized Batching**

```python
# Calculate optimal batch size based on available memory
available_memory_gb = 32  # GB
memory_per_sample_gb = 0.1  # GB per sample (adjust based on your data)

optimal_batch_size = int((available_memory_gb * 0.8) / memory_per_sample_gb)  # Leave 20% buffer

processor = BatchProcessor(
    batch_size=min(optimal_batch_size, 100),  # Cap at 100 samples per batch
    n_jobs=min(16, multiprocessing.cpu_count() - 2),  # Leave cores for OS
)
```

### Performance Optimization Tips

#### 1. **Sample Grouping Strategy**

```python
# Group samples by chromosome to optimize I/O
def group_samples_by_chromosome(samples, variants):
    chr_groups = {}
    for sample in samples:
        # Analyze which chromosomes each sample's variants are on
        sample_chrs = set(variant.chrom for variant in variants)
        primary_chr = max(sample_chrs, key=lambda x: sum(1 for v in variants if v.chrom == x))

        if primary_chr not in chr_groups:
            chr_groups[primary_chr] = []
        chr_groups[primary_chr].append(sample)

    return chr_groups

# Process each chromosome group separately
for chr_group in chromosome_groups.values():
    process_sample_group(chr_group, variants_for_chromosome)
```

#### 2. **Variant Pre-filtering**

```bash
# Pre-filter variants to reduce processing load
bcftools view large_panel.vcf -R target_regions.bed > filtered_panel.vcf

# Process only relevant variants
gbcms count run \
    --bam-fof samples.txt \
    --vcf filtered_panel.vcf \
    --thread 32
```

#### 3. **Progressive Processing**

```python
# Process in phases for very large sample sets
def process_in_phases(samples, variants, samples_per_phase=100):
    results = []
    for i in range(0, len(samples), samples_per_phase):
        phase_samples = samples[i:i+samples_per_phase]

        # Process this phase
        phase_results = process_sample_batch(phase_samples, variants)

        # Save intermediate results
        save_intermediate_results(phase_results, f"phase_{i//samples_per_phase}.txt")

        results.extend(phase_results)

        # Optional: memory cleanup
        gc.collect()

    return results
```

### Troubleshooting High Sample Counts

#### 1. **Memory Issues**

**Problem**: Out of memory with 500+ samples
**Solutions**:
```python
# Use smaller batch sizes
processor = BatchProcessor(batch_size=25, n_jobs=8)

# Process samples in phases
results = process_in_phases(all_samples, variants, samples_per_phase=200)

# Monitor memory usage
import psutil
memory_usage = psutil.virtual_memory().percent
if memory_usage > 80:
    # Reduce batch size or switch to Ray
```

#### 2. **I/O Bottlenecks**

**Problem**: Slow with many BAM files
**Solutions**:
```bash
# Use threading backend for I/O-heavy workloads
gbcms count run --backend threading --thread 8

# Pre-load BAM indices
for bam in all_bam_files:
    ensure_bam_indexed(bam)
```

#### 3. **Network Issues (Ray)**

**Problem**: Ray workers timing out
**Solutions**:
```python
# Increase Ray timeouts for large clusters
ray.init(address=cluster_address, timeout=300)  # 5 minute timeout

# Use smaller worker counts initially
processor = ParallelProcessor(n_jobs=32, backend='ray')  # Start smaller
```

### Real-World Example: 500 Clinical Samples

```bash
# Scenario: Process 500 clinical samples against a 500-gene panel

# Step 1: Prepare sample list
cat > clinical_samples.txt << EOF
patient001	patient001.bam
patient002	patient002.bam
# ... 498 more samples
patient500	patient500.bam
EOF

# Step 2: Process with optimized settings
gbcms count run \
    --fasta GRCh38.fa \
    --bam-fof clinical_samples.txt \
    --vcf clinical_panel.vcf \
    --output clinical_results.txt \
    --thread 32 \
    --backend joblib \
    --max-block-size 1000 \
    --positive-count \
    --fragment-count

# Expected runtime: 20-40 minutes
# Output: ~2.5GB text file with all sample counts
```

### Summary for High Sample Counts

| 100-200 | joblib | 16-24 threads | 5-15 min | 4-8 GB |
| 200-500 | joblib | 24-32 threads | 15-45 min | 8-16 GB |
| 500-1000 | Ray | 48-96 workers | 30-90 min | Distributed |
| 1000+ | Ray | 96-200+ workers | 1-3 hours | Distributed |

**Key Takeaway**: For most clinical/research scenarios with 100s-1000s of samples, **joblib with 16-32 threads** provides the best balance of performance, simplicity, and reliability.

## 📊 Decision Tree

```
Need parallelization?
├── Small dataset (<10K variants)
│   └── Use single-threaded or joblib with 2-4 threads
│
├── Medium dataset (10K-100K variants)
│   └── Use joblib with 8-16 threads
│       ├── I/O heavy? → threading backend
│       └── CPU heavy? → loky backend
│
├── Large dataset (100K-1M variants)
│   └── Use joblib with 16-32 threads
│       ├── Many samples? → Consider batching
│       └── Memory issues? → Use BatchProcessor
│
└── Massive dataset (>1M variants)
    └── Use Ray for distributed processing
        ├── Single machine? → Ray with multiple workers
        └── Multi-node cluster? → Full Ray cluster
```

**For High Sample Counts (100s-1000s of samples):**
- **100-500 samples**: joblib with 16-32 threads (15-45 minutes)
- **500-1000+ samples**: Ray with 48-96+ workers (30 minutes - 2 hours)

## ⚡ Performance Tips

### joblib Optimization

1. **Choose right backend:**
{{ ... }}
   # CPU-intensive (counting): loky
   # I/O-intensive (loading): threading
   # Maximum parallelism: multiprocessing
   ```

2. **Batch size matters:**
   ```python
   # For large datasets: 1000-5000 variants per batch
   # For small datasets: 100-500 variants per batch
   batch_size = min(1000, len(variants) // n_jobs)
   ```

3. **Memory management:**
   ```python
   # Process in chunks to avoid OOM
   chunk_size = 10000
   for i in range(0, len(variants), chunk_size):
       chunk = variants[i:i+chunk_size]
       results.extend(parallel_map(func, chunk, n_jobs=8))
   ```

### Ray Optimization

1. **Worker configuration:**
   ```python
   # One actor per physical core for best performance
   n_actors = min(64, total_cluster_cores)
   ```

2. **Memory per worker:**
   ```python
   # Ray reserves memory per worker
   # Leave 2-4GB per worker for OS + overhead
   memory_per_worker = (total_memory - 8) // n_workers
   ```

3. **Fault tolerance:**
   ```python
   # Ray automatically handles worker failures
   # Results are collected even if some workers fail
   # Failed tasks are retried automatically
   ```

## 🔧 Troubleshooting

### joblib Issues

**Problem**: High memory usage
**Solution**: Use BatchProcessor with smaller batch sizes

**Problem**: Poor performance with I/O
**Solution**: Use threading backend instead of loky

**Problem**: GIL limitations
**Solution**: Use multiprocessing backend for CPU-only tasks

### Ray Issues

**Problem**: Workers failing
**Solution**: Check Ray dashboard for worker status and logs

**Problem**: Network timeouts
**Solution**: Increase Ray timeout settings for large clusters

**Problem**: Memory pressure
**Solution**: Reduce workers per node or increase node memory

## 📚 Additional Resources

- **[Advanced Features Guide](ADVANCED_FEATURES.md)** - Detailed API documentation
- **[Architecture Overview](ARCHITECTURE.md)** - How parallelization fits into gbcms
- **[Performance Benchmarks](ADVANCED_FEATURES.md#performance-benchmarks)** - Real-world performance data
- **[Ray Documentation](https://docs.ray.io/)** - Official Ray documentation
- **[joblib Documentation](https://joblib.readthedocs.io/)** - Official joblib documentation

---

**Need help choosing?** Start with joblib - it's simpler and works for 90% of use cases. Only use Ray when you need massive scale or fault tolerance.
