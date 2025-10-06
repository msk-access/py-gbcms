# Parallelization Guide: Multi-Level Optimization Strategy

Complete guide to understanding and optimizing gbcms parallel processing for your genomics workloads.

## 🏗️ **Understanding gbcms Parallel Processing**

gbcms uses a sophisticated **multi-level parallel processing strategy** that automatically adapts to your input:

### **Three-Level Architecture**

| **Level** | **Granularity** | **Parallelization** | **Primary Benefit** |
|-----------|----------------|--------------------|-------------------|
| **🔬 BAM-Level** | **Inter-Sample** | Multiple BAMs → Multiple threads | Scales with sample count |
| **📦 Block-Level** | **Intra-Sample** | Variant blocks → Thread pools | Optimizes large variant lists |
| **⚡ Algorithm-Level** | **Computational** | Numba JIT + vectorization | Maximizes counting performance |

### **How It Works**

**🎯 **Multiple BAM Scenario** (Multiple Sample Analysiss)**:


**🧬 **Single BAM Scenario** (Large Variant Lists)**:


---

## 🎯 **Scenario-Based Recommendations**

### **🏥 **Multiple Sample Processing**

**Best For**: Sample collections, large-scale studies, genomics analysis


**💡 **Why This Works**:
- **BAM-Level**: Parallelizes across samples (if multiple BAMs)
- **Block-Level**: Handles large variant lists efficiently
- **Threading**: Optimal for I/O-bound genomics workloads
- **Expected Speedup**: 4-16× depending on sample count

### **🧬 **Single Sample Analysis (Single Large BAM)**

**Best For**: Deep sequencing, large variant calling studies


**💡 **Why This Works**:
- **Block-Level**: Optimizes processing of 10K+ variants
- **High Thread Count**: Maximizes intra-sample parallelization
- **Fragment Weighting**: Better for deep sequencing data
- **Expected Speedup**: 8-32× for large variant lists

### **💻 **Development & Testing**

**Best For**: Algorithm development, small datasets, debugging


**💡 **Why This Works**:
- **Threading Backend**: Fast startup, shared memory
- **Smaller Blocks**: Quick feedback for development
- **Lower Thread Count**: Laptop-friendly resource usage

### **🏭 **High-Throughput Production**

**Best For**: Genomics laboratories, core facilities, production pipelines


**💡 **Why This Works**:
- **High Thread Count**: Maximizes cluster utilization
- **Large Blocks**: Minimizes I/O overhead
- **Strict Filtering**: Quality control for production
- **Expected Speedup**: 16-64× for large batches

---

## 📊 **Performance Expectations**

### **Scaling Characteristics**

| **Input Size** | **Recommended Threads** | **Expected Speedup** | **Memory Usage** |
|---------------|------------------------|---------------------|------------------|
| **1-4 BAMs** | **4-8 threads** | **3-6×** | **Low-Moderate** |
| **5-20 BAMs** | **8-16 threads** | **6-12×** | **Moderate** |
| **20+ BAMs** | **16-32 threads** | **12-24×** | **High** |
| **1 BAM, 10K+ variants** | **16-32 threads** | **8-24×** | **Moderate-High** |

### **Hardware Recommendations**

**💻 **Workstation/Laptop** (8-16 cores)**:


**🖥️ **Small Server** (16-32 cores)**:


**🏢 **HPC Cluster** (32-128 cores)**:


---

## ⚙️ **Configuration Options Explained**

### **Thread Configuration ()**

| **Setting** | **Use Case** | **Performance** | **Resource Usage** |
|-------------|-------------|----------------|-------------------|
| **1-4** | Development, small datasets | **Baseline** | **Minimal** |
| **4-8** | Workstations, single samples | **Good** | **Low** |
| **8-16** | Small servers, moderate batches | **Excellent** | **Moderate** |
| **16-32** | HPC nodes, large batches | **Outstanding** | **High** |
| **32+** | Large clusters, massive datasets | **Maximum** | **Very High** |

### **Backend Selection ()**

**🔧 **joblib (Recommended)**:


**🧵 **threading**:


**⚡ **multiprocessing**:


### **Block Size Configuration**

**📦 **Max Block Size ()**:
- **Small (100-500)**: More blocks, better load balancing, more I/O
- **Medium (500-2000)**: Balanced I/O and memory usage
- **Large (2000+)**: Fewer blocks, less I/O, more memory per block

**📏 **Max Block Distance ()**:
- **Small (10K-50K)**: More blocks, better for clustered variants
- **Medium (50K-200K)**: Good balance for most genomes
- **Large (200K+)**: Fewer blocks, better for sparse variants

---

## 🚀 **Quick Start Examples**

### **Basic Usage**


### **Multiple Sample Analysis**


### **Single Sample Analysis**


### **Development**


---

## 🔧 **Troubleshooting Performance**

### **Common Issues & Solutions**

**🐌 **Slow Performance**:


**💾 **High Memory Usage**:


**⚠️ **Errors or Crashes**:


### **Performance Monitoring**

**📊 **Monitor Resource Usage**:


**🎯 **Expected Performance**:
- **Small Dataset (1-5 samples)**: 2-8× speedup
- **Medium Dataset (5-20 samples)**: 8-16× speedup
- **Large Dataset (20+ samples)**: 16-32× speedup
- **Memory Usage**: 2-8GB typical, 8-32GB for large datasets

---

## 💡 **Best Practices Summary**

1. **🎯 **Start Simple**: Use  for most cases
2. **📈 **Scale Threads**: Match thread count to available CPU cores
3. **⚖️ **Balance Load**: Use  for balanced I/O/memory
4. **🔍 **Monitor**: Watch CPU, memory, and I/O during processing
5. **⚡ **Optimize**: Adjust based on your specific data characteristics

**Remember**: gbcms automatically adapts its parallelization strategy based on your input - whether you have multiple BAMs (BAM-level parallelization) or single large BAMs (block-level parallelization). The multi-level approach ensures optimal performance across different genomics workloads! 🚀
