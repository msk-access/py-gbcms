# Variant Counting

How py-gbcms counts reads for each variant type — with visual examples.

## Overview

For every variant, py-gbcms fetches reads overlapping the variant position from each BAM file, applies [read filters](#read-filters), and then dispatches to a **type-specific** allele checker. Each read is classified as supporting the **reference** allele, the **alternate** allele, or **neither**.

```mermaid
flowchart LR
    Fetch["📥 Fetch reads\nat variant locus"] --> Filter["🔍 Apply\nread filters"]
    Filter --> Dispatch{"Variant\ntype?"}
    Dispatch -->|SNP| SNP["check_snp"]
    Dispatch -->|Insertion| INS["check_insertion"]
    Dispatch -->|Deletion| DEL["check_deletion"]
    Dispatch -->|MNP| MNP["check_mnp"]
    Dispatch -->|Complex| CPX["check_complex"]
    SNP --> Count["📊 Update\ncounts"]
    INS --> Count
    DEL --> Count
    MNP --> Count
    CPX --> Count

    style Fetch fill:#4a90d9,color:#fff
    style Filter fill:#e67e22,color:#fff
    style Count fill:#27ae60,color:#fff
```

---

## Read Filter Pipeline

Before any allele checking begins, every read passes through a **filter cascade**. Reads that fail any enabled filter are discarded. The order matches the Rust engine implementation.

```mermaid
flowchart TD
    Read["📖 Read from BAM"] --> F1{"🔴 Duplicate?\n--filter-duplicates"}
    F1 -->|"Yes (default: ON)"| Drop["❌ Discard read"]
    F1 -->|No| F2{"Secondary?\n--filter-secondary"}
    F2 -->|"Yes & filter ON"| Drop
    F2 -->|No or filter OFF| F3{"Supplementary?\n--filter-supplementary"}
    F3 -->|"Yes & filter ON"| Drop
    F3 -->|No or filter OFF| F4{"QC Failed?\n--filter-qc-failed"}
    F4 -->|"Yes & filter ON"| Drop
    F4 -->|No or filter OFF| F5{"Improper pair?\n--filter-improper-pair"}
    F5 -->|"Yes & filter ON"| Drop
    F5 -->|No or filter OFF| F6{"Contains indel?\n--filter-indel"}
    F6 -->|"Yes & filter ON"| Drop
    F6 -->|No or filter OFF| F7{"MAPQ ≥ 20?\n--min-mapq"}
    F7 -->|"No (below threshold)"| Drop
    F7 -->|"Yes"| Pass["✅ Pass to\nallele checker"]

    style Read fill:#3498db,color:#fff
    style Drop fill:#e74c3c,color:#fff
    style Pass fill:#27ae60,color:#fff
```

| Filter | CLI Flag | Default | SAM Flag |
|:-------|:---------|:--------|:---------|
| Duplicates | `--filter-duplicates` | **On** | `0x400` |
| Secondary | `--filter-secondary` | Off | `0x100` |
| Supplementary | `--filter-supplementary` | Off | `0x800` |
| QC Failed | `--filter-qc-failed` | Off | `0x200` |
| Improper Pair | `--filter-improper-pair` | Off | `0x2` (inverted) |
| Indel reads | `--filter-indel` | Off | CIGAR-based |
| MAPQ threshold | `--min-mapq` | **20** | — |
| BASEQ threshold | `--min-baseq` | **0** | — (per-type) |

!!! note "Filter Non-Primary"
    The original GBCMS has a single `--filter_non_primary` flag. py-gbcms splits this into `--filter-secondary` and `--filter-supplementary` for finer control. Both default to **off**, matching the original behavior.

---

## Variant Types

### SNP (Single Nucleotide Polymorphism)

A single base substitution — the simplest and most common variant type.

| Property | Value |
|:---------|:------|
| Detection | `len(REF) == 1 && len(ALT) == 1` |
| Position | 0-based index of the substituted base |
| Quality check | Base quality at the position must meet `--min-baseq` |

#### Algorithm Flow

```mermaid
flowchart TD
    Start["🧬 SNP Check\npos, REF, ALT"] --> Walk["Walk CIGAR to find\nread position at\nvariant.pos"]
    Walk --> Found{"Position\nfound?"}
    Found -->|No| Neither["Neither\n(read doesn't cover)"]
    Found -->|Yes| BQ{"Base quality\n≥ min_baseq?"}
    BQ -->|No| Neither2["Neither\n(low quality)"]
    BQ -->|Yes| Compare["Compare base to\nREF and ALT\n(case-insensitive)"]
    Compare --> IsRef{"base == REF?"}
    IsRef -->|Yes| Ref["✅ REF"]
    IsRef -->|No| IsAlt{"base == ALT?"}
    IsAlt -->|Yes| Alt["🔴 ALT"]
    IsAlt -->|No| Neither3["Neither\n(third allele)"]

    style Start fill:#9b59b6,color:#fff
    style Ref fill:#27ae60,color:#fff
    style Alt fill:#e74c3c,color:#fff
    style Neither fill:#95a5a6,color:#fff
    style Neither2 fill:#95a5a6,color:#fff
    style Neither3 fill:#95a5a6,color:#fff
```

#### Visual Example

```
Variant: chr1:100 C→T (SNP)

Reference: 5'─ ...G  A  T  C  G  A  T  C  G  A... ─3'
                          98 99 100 101
                               ▲
                           variant pos

Read 1:    5'─ ...G  A  T [T] G  A  T  C... ─3'   → ALT ✅
                               ↑
                          base=T matches ALT

Read 2:    5'─ ...G  A  T [C] G  A  T  C... ─3'   → REF ✅
                               ↑
                          base=C matches REF

Read 3:    5'─ ...G  A  T [A] G  A  T  C... ─3'   → Neither
                               ↑
                          base=A ≠ REF or ALT
```

---

### Insertion

Bases inserted after an **anchor** position. The anchor is the last reference base before the inserted sequence.

| Property | Value |
|:---------|:------|
| Detection | `len(REF) == 1 && len(ALT) > 1` |
| Position | 0-based index of the **anchor** base |
| Quality check | None (CIGAR-structural check only) |

#### Algorithm Flow

```mermaid
flowchart TD
    Start["🧬 Insertion Check\nanchor_pos, REF=A, ALT=ATG"] --> Walk["Walk CIGAR ops\nto find match block\ncontaining anchor"]
    Walk --> Found{"Anchor\nfound in\nmatch block?"}
    Found -->|No| Neither["Neither"]
    Found -->|Yes| AtEnd{"Anchor is\nlast base of\nmatch block?"}
    AtEnd -->|No| Ref["✅ REF\n(anchor in middle\nof match = no insertion)"]
    AtEnd -->|Yes| NextOp{"Next CIGAR\nop is Ins?"}
    NextOp -->|No| Ref2["✅ REF\n(no insertion follows)"]
    NextOp -->|Yes| LenCheck{"Insertion length\n== len(ALT) - 1?"}
    LenCheck -->|No| Neither2["Neither\n(wrong length)"]
    LenCheck -->|Yes| SeqCheck{"Inserted sequence\n== ALT[1..]?"}
    SeqCheck -->|No| Neither3["Neither\n(wrong sequence)"]
    SeqCheck -->|Yes| Alt["🔴 ALT"]

    style Start fill:#9b59b6,color:#fff
    style Ref fill:#27ae60,color:#fff
    style Ref2 fill:#27ae60,color:#fff
    style Alt fill:#e74c3c,color:#fff
    style Neither fill:#95a5a6,color:#fff
    style Neither2 fill:#95a5a6,color:#fff
    style Neither3 fill:#95a5a6,color:#fff
```

#### Visual Example

```
Variant: chr1:100 A→ATG (insertion of TG after anchor A)

Reference:     5'─ ...C  G  A ── C  G  T  A... ─3'
                          99 100  101 102
                               ▲
                          anchor pos

Read 1 (ALT):  CIGAR = 5M 2I 5M
               5'─ ...C  G  A [T  G] C  G  T  A... ─3'
                               └──┘
                          inserted bases match ALT → ALT ✅

Read 2 (REF):  CIGAR = 10M
               5'─ ...C  G  A  C  G  T  A... ─3'
                               ↑
                          no insertion after anchor → REF ✅

Read 3 (other): CIGAR = 5M 1I 5M
               5'─ ...C  G  A [C] C  G  T  A... ─3'
                               ↑
                          1bp insertion ≠ expected 2bp → Neither
```

!!! tip "Why anchor-based?"
    In VCF format, insertions are represented as `REF=A, ALT=ATG` where the first base `A` is the anchor — it's not part of the inserted sequence. py-gbcms uses this anchor to locate where the insertion should occur in the CIGAR string.

---

### Deletion

Bases deleted after an **anchor** position. The logic mirrors insertion but looks for `Del` CIGAR operations.

| Property | Value |
|:---------|:------|
| Detection | `len(REF) > 1 && len(ALT) == 1` |
| Position | 0-based index of the **anchor** base |
| Quality check | None (CIGAR-structural check only) |

#### Algorithm Flow

```mermaid
flowchart TD
    Start["🧬 Deletion Check\nanchor_pos, REF=ACG, ALT=A"] --> Walk["Walk CIGAR ops\nto find match block\ncontaining anchor"]
    Walk --> Found{"Anchor\nfound in\nmatch block?"}
    Found -->|No| Neither["Neither"]
    Found -->|Yes| AtEnd{"Anchor is\nlast base of\nmatch block?"}
    AtEnd -->|No| Ref["✅ REF\n(anchor in middle\nof match = no deletion)"]
    AtEnd -->|Yes| NextOp{"Next CIGAR\nop is Del?"}
    NextOp -->|No| Ref2["✅ REF\n(no deletion follows)"]
    NextOp -->|Yes| LenCheck{"Deletion length\n== len(REF) - 1?"}
    LenCheck -->|No| Neither2["Neither\n(wrong length)"]
    LenCheck -->|Yes| Alt["🔴 ALT"]

    style Start fill:#9b59b6,color:#fff
    style Ref fill:#27ae60,color:#fff
    style Ref2 fill:#27ae60,color:#fff
    style Alt fill:#e74c3c,color:#fff
    style Neither fill:#95a5a6,color:#fff
    style Neither2 fill:#95a5a6,color:#fff
```

#### Visual Example

```
Variant: chr1:100 ACG→A (deletion of CG after anchor A)

Reference:     5'─ ...T  G  A  C  G  T  A  C... ─3'
                          99 100 101 102
                               ▲
                          anchor pos

Read 1 (ALT):  CIGAR = 5M 2D 5M
               5'─ ...T  G  A  ──  ──  T  A  C... ─3'
                               └─────┘
                          2bp deletion matches len(REF)-1 → ALT ✅

Read 2 (REF):  CIGAR = 12M
               5'─ ...T  G  A  C  G  T  A  C... ─3'
                               ↑
                          no deletion after anchor → REF ✅

Read 3 (other): CIGAR = 5M 3D 5M
               5'─ ...T  G  A  ──  ──  ──  A  C... ─3'
                               └────────┘
                          3bp deletion ≠ expected 2bp → Neither
```

---

### MNP (Multi-Nucleotide Polymorphism)

Multiple adjacent bases substituted simultaneously. Think of it as multiple SNPs happening at consecutive positions.

| Property | Value |
|:---------|:------|
| Detection | `len(REF) == len(ALT) && len(REF) > 1` |
| Position | 0-based index of the first substituted base |
| Quality check | **Every** base in the MNP region must meet `--min-baseq` |

#### Algorithm Flow

```mermaid
flowchart TD
    Start["🧬 MNP Check\npos, REF=AT, ALT=CG"] --> Find["Find read position\nof first base"]
    Find --> Found{"Position\nfound?"}
    Found -->|No| Neither["Neither"]
    Found -->|Yes| Cover{"Read covers\nentire MNP\nregion?"}
    Cover -->|No| Neither2["Neither"]
    Cover -->|Yes| Loop["For each position\nin MNP region:"]
    Loop --> BQ{"Base quality\n≥ min_baseq?"}
    BQ -->|No| Neither3["Neither\n(any low-quality base\nfails entire MNP)"]
    BQ -->|Yes| Track["Compare to REF[i]\nand ALT[i]"]
    Track --> More{"More\npositions?"}
    More -->|Yes| Loop
    More -->|No| Contig{"Contiguity check:\nno indels within\nMNP region?"}
    Contig -->|"Indel found"| Neither4["Neither"]
    Contig -->|"Contiguous"| Final{"All bases\nmatch?"}
    Final -->|"All match ALT"| Alt["🔴 ALT"]
    Final -->|"All match REF"| Ref["✅ REF"]
    Final -->|"Mixed"| Neither5["Neither\n(partial match)"]

    style Start fill:#9b59b6,color:#fff
    style Ref fill:#27ae60,color:#fff
    style Alt fill:#e74c3c,color:#fff
    style Neither fill:#95a5a6,color:#fff
    style Neither2 fill:#95a5a6,color:#fff
    style Neither3 fill:#95a5a6,color:#fff
    style Neither4 fill:#95a5a6,color:#fff
    style Neither5 fill:#95a5a6,color:#fff
```

#### Visual Example

```
Variant: chr1:100 AT→CG (2-base MNP)

Reference:     5'─ ...G  A  T  A  T  G  C... ─3'
                          99 100 101
                              ▲───▲
                           MNP region

Read 1 (ALT):  5'─ ...G  A  T [C  G] G  C... ─3'   → ALT ✅
                              └──┘
                         both bases match ALT: C=C ✓, G=G ✓

Read 2 (REF):  5'─ ...G  A  T [A  T] G  C... ─3'   → REF ✅
                              └──┘
                         both bases match REF: A=A ✓, T=T ✓

Read 3 (mixed): 5'─ ...G  A  T [C  T] G  C... ─3'  → Neither
                              └──┘
                         C=ALT[0] ✓, but T=REF[1] ✗  (partial match)

Read 4 (indel): CIGAR = 3M 1I 3M
                5'─ ...G  A  T [C] X [G] G  C... ─3' → Neither
                                   ↑
                         indel within MNP region → fails contiguity check
```

!!! warning "Strict matching"
    MNP matching is **all-or-nothing**: every base must match either REF or ALT. A read with `C T` at a `AT→CG` variant (matching the first ALT base but second REF base) is classified as **neither** — it's not counted for either allele.

---

### Complex (Indel + Substitution)

Variants where REF and ALT differ in both sequence **and** length. This is the catch-all category that uses a sophisticated **haplotype reconstruction** algorithm.

| Property | Value |
|:---------|:------|
| Detection | Fallback for all other combinations |
| Position | 0-based index of the first reference base |
| Quality check | Minimum base quality across all reconstructed bases must meet `--min-baseq` |

#### Algorithm Flow — Haplotype Reconstruction

```mermaid
flowchart TD
    Start["🧬 Complex Check\npos=10, REF=ACG, ALT=TT"] --> Region["Define genomic region\n[pos, pos + len(REF))\n= [10, 13)"]
    Region --> Init["Initialize:\n• reconstructed_seq = []\n• min_qual = 255"]
    Init --> Walk["Walk each CIGAR op"]
    Walk --> OpType{"CIGAR op\ntype?"}

    OpType -->|"M / = / X\n(Match)"| Match["Does op overlap\nwith [10, 13)?"]
    Match -->|Yes| AppendM["Append overlapping\nbases from read.\nTrack min quality."]
    Match -->|No| Next["Next CIGAR op"]
    AppendM --> Next

    OpType -->|"I (Insertion)"| InsCheck["Is ref_pos within\n[10, 13]?"]
    InsCheck -->|Yes| AppendI["Append inserted\nbases from read"]
    InsCheck -->|No| Next
    AppendI --> Next

    OpType -->|"D / N\n(Del/Skip)"| AdvRef["Advance ref_pos\nonly (no bases\nappended)"]
    AdvRef --> Next

    OpType -->|"S (Soft clip)"| AdvRead["Advance read_pos\nonly"]
    AdvRead --> Next

    OpType -->|"H / P\n(Hard/Pad)"| Next

    Next --> More{"More\nCIGAR ops?"}
    More -->|Yes| Walk
    More -->|No| Compare["Compare\nreconstructed\nsequence"]

    Compare --> CmpAlt{"reconstructed\n== ALT?"}
    CmpAlt -->|Yes| QualA{"min_qual\n≥ min_baseq?"}
    QualA -->|Yes| Alt["🔴 ALT"]
    QualA -->|No| Neither["Neither\n(low quality)"]

    CmpAlt -->|No| CmpRef{"reconstructed\n== REF?"}
    CmpRef -->|Yes| QualR{"min_qual\n≥ min_baseq?"}
    QualR -->|Yes| Ref["✅ REF"]
    QualR -->|No| Neither2["Neither\n(low quality)"]
    CmpRef -->|No| Neither3["Neither\n(no match)"]

    style Start fill:#9b59b6,color:#fff
    style Region fill:#3498db,color:#fff
    style Ref fill:#27ae60,color:#fff
    style Alt fill:#e74c3c,color:#fff
    style Neither fill:#95a5a6,color:#fff
    style Neither2 fill:#95a5a6,color:#fff
    style Neither3 fill:#95a5a6,color:#fff
```

#### Worked Example: Step-by-Step Reconstruction

Here's a concrete example showing how the engine reconstructs a read's sequence for a complex variant.

```
Variant: chr1:10 ACG→TT (complex: 3bp REF → 2bp ALT)
Region:  [10, 13)   (0-based, half-open)
```

**Read supporting ALT** — CIGAR: `10M 1D 1I 8M`

```
Step-by-step CIGAR walk:

  CIGAR op    ref_pos  read_pos  Overlap with [10,13)?  Action
  ─────────   ───────  ────────  ─────────────────────  ──────────────────────
  10M         0→10     0→10      [10,10) = none         (just advance)
  ⚠️ Wait — let's reconsider. 10M covers ref [0,10), read [0,10).
     overlap with [10,13) = none yet.

  Actually, the 10M goes from ref 0 to 10 (exclusive).
  After 10M: ref_pos=10, read_pos=10.

  1D          10→11    (no read)  ref_pos 10 is in [10,13)
                                  But Del → no bases appended, ref advances
              After: ref_pos=11, read_pos=10

  1I          11       10→11     ref_pos=11 is in [10,13]
                                  Append 1 inserted base from read[10]
                                  reconstructed = "T"
              After: ref_pos=11, read_pos=11

  8M          11→19    11→19     overlap [11,13) = 2 bases
                                  Append read[11] and read[12]
                                  reconstructed = "T" + "T" = "TT"
              After: ref_pos=19, read_pos=19

  ─────────────────────────────────────────────────────────────────
  Result: reconstructed = "TT"
          Compare to ALT = "TT" → ✅ Match! → ALT
```

**Read supporting REF** — CIGAR: `20M`

```
  CIGAR op    ref_pos  read_pos  Overlap with [10,13)?  Action
  ─────────   ───────  ────────  ─────────────────────  ──────────────
  20M         0→20     0→20      [10,13) = 3 bases
                                  Append read[10], read[11], read[12]
                                  reconstructed = "ACG"

  ─────────────────────────────────────────────────────────────────
  Result: reconstructed = "ACG"
          Compare to ALT = "TT"  → length 3 ≠ 2 → no match
          Compare to REF = "ACG" → ✅ Match! → REF
```

!!! info "Auto-detection"
    If a variant's type string is unknown, py-gbcms auto-detects MNPs (equal-length REF/ALT with length > 1) and falls back to the complex haplotype reconstruction for everything else.

---

## Counting Levels

After allele classification, py-gbcms tracks counts at **two levels**: individual reads and collapsed fragments.

### Read-Level vs Fragment-Level

```mermaid
flowchart TD
    subgraph ReadLevel["📖 Read-Level Counting"]
        R1["Read 1 (fwd) → ALT"]
        R2["Read 2 (rev) → ALT"]
        R3["Read 3 (fwd) → REF"]
        R4["Read 4 (rev) → REF"]
        R5["Read 5 (fwd) → ALT"]
        R6["Read 6 (rev) → REF"]
    end

    subgraph FragLevel["🧬 Fragment-Level Counting"]
        direction TB
        F1["Fragment A\nR1+R2 → ALT\norientation: fwd (R1)"]
        F2["Fragment B\nR3+R4 → REF\norientation: fwd (R3)"]
        F3["Fragment C\nR5 only → ALT\norientation: fwd"]
        F4["Fragment D\nR6 only → REF\norientation: rev"]
    end

    R1 --> F1
    R2 --> F1
    R3 --> F2
    R4 --> F2
    R5 --> F3
    R6 --> F4

    subgraph Results["📊 Final Counts"]
        ReadCounts["Read Level:\nDP=6, RD=3, AD=3\nDP_fwd=3, DP_rev=3"]
        FragCounts["Fragment Level:\nDPF=4, RDF=2, ADF=2\nRDF_fwd=1, RDF_rev=1"]
    end

    F1 --> FragCounts
    F2 --> FragCounts
    F3 --> FragCounts
    F4 --> FragCounts

    style ReadLevel fill:#3498db15,stroke:#3498db
    style FragLevel fill:#9b59b615,stroke:#9b59b6
    style Results fill:#27ae6015,stroke:#27ae60
```

### Read-Level Metrics

Each read is counted independently.

| Metric | Description |
|:-------|:------------|
| **DP** | Total depth (reads supporting REF or ALT) |
| **RD** / **AD** | Reference / Alternate read counts |
| **DP_fwd** / **DP_rev** | Strand-specific total depth |
| **RD_fwd** / **RD_rev** | Strand-specific reference counts |
| **AD_fwd** / **AD_rev** | Strand-specific alternate counts |

### Fragment Counting Algorithm

Fragment counting collapses read pairs into a single observation per fragment.

```mermaid
flowchart TD
    Start["For each read passing filters:"] --> Track["Track by fragment name (QNAME):\n• Which allele(s) each read supports\n• Orientation of each read"]
    Track --> Done["After all reads processed:"]
    Done --> ForEach["For each unique fragment:"]
    ForEach --> HasBoth{"Both reads\npresent?"}
    HasBoth -->|Yes| Orient1["Orientation = Read 1's direction\n(forward or reverse)"]
    HasBoth -->|No| Orient2["Orientation = single read's direction"]
    Orient1 --> Allele
    Orient2 --> Allele
    Allele["Check alleles:\n• has_ref? → increment RDF\n• has_alt? → increment ADF\n• Always increment DPF"]
    Allele --> Strand["Update strand-specific\ncounts based on\nfragment orientation"]

    style Start fill:#3498db,color:#fff
```

| Metric | Description |
|:-------|:------------|
| **DPF** | Fragment depth |
| **RDF** / **ADF** | Reference / Alternate fragment counts |
| **RDF_fwd** / **RDF_rev** | Strand-specific reference fragment counts |
| **ADF_fwd** / **ADF_rev** | Strand-specific alternate fragment counts |

!!! tip "Why fragment counting matters for cfDNA"
    In cfDNA sequencing (like MSK-ACCESS), the same DNA fragment is often sequenced from both ends. Counting reads would double-count each fragment. Fragment-level counting gives a more accurate picture of the number of **unique molecules** supporting each allele.

---

### Strand Bias (Fisher's Exact Test)

Strand bias detects when an allele is disproportionately supported by reads on one strand — a common sequencing artifact.

```mermaid
flowchart LR
    subgraph Table["2×2 Contingency Table"]
        direction TB
        Header["| · | Forward | Reverse |"]
        RefRow["| REF | a (RD_fwd) | b (RD_rev) |"]
        AltRow["| ALT | c (AD_fwd) | d (AD_rev) |"]
    end

    Table --> Fisher["Fisher's Exact\nTest"]
    Fisher --> Pval["SB_pval\n(p-value)"]
    Fisher --> OR["SB_OR\n(odds ratio)"]

    Pval --> Interpret{"p < 0.05?"}
    Interpret -->|Yes| Bias["⚠️ Possible\nstrand bias"]
    Interpret -->|No| NoBias["✅ No evidence\nof strand bias"]

    style Bias fill:#e74c3c,color:#fff
    style NoBias fill:#27ae60,color:#fff
```

Computed at **both** levels:

| Metric | Level | Description |
|:-------|:------|:------------|
| **SB_pval** / **SB_OR** | Read | Strand bias from individual reads |
| **FSB_pval** / **FSB_OR** | Fragment | Strand bias from collapsed fragments |

!!! example "Strand bias example"
    If a variant has `AD_fwd=15, AD_rev=1`, that's suspicious — almost all ALT-supporting reads are on the forward strand. Fisher's test would yield a low p-value, flagging this as a potential artifact rather than a true variant.

---

## Full Pipeline: End-to-End Example

Here's how a single variant is processed through the complete pipeline:

```mermaid
sequenceDiagram
    participant CLI as CLI (Python)
    participant Pipeline as Pipeline
    participant Reader as VCF/MAF Reader
    participant Kernel as Coordinate Kernel
    participant Rust as Rust Engine
    participant BAM as BAM File

    CLI->>Pipeline: run(config)
    Pipeline->>Reader: load variants
    Reader->>Kernel: normalize coordinates
    Kernel-->>Reader: 0-based Variant objects
    Reader-->>Pipeline: List[Variant]

    loop For each BAM sample
        Pipeline->>Rust: count_bam(bam, variants, filters)
        loop For each variant (parallel via Rayon)
            Rust->>BAM: fetch(chrom, pos, pos+1)
            BAM-->>Rust: Iterator of reads
            loop For each read
                Note over Rust: Apply filter cascade
                Note over Rust: Dispatch to type checker
                Note over Rust: Update read + fragment counts
            end
            Note over Rust: Compute Fisher strand bias
        end
        Rust-->>Pipeline: Vec[BaseCounts]
    end

    Pipeline->>Pipeline: Write output (VCF/MAF)
```

---

## Comparison with Original GBCMS

| Feature | Original GBCMS | py-gbcms |
|:--------|:---------------|:---------|
| Counting algorithm | Region-based chunking, position matching | Per-variant CIGAR traversal |
| Complex variants | Optional via `--generic_counting` | Always uses haplotype reconstruction |
| MNP handling | Not explicit | Dedicated `check_mnp` with contiguity check |
| Fragment counting | Optional (`--fragment_count`) | Always computed |
| Positive strand counts | Optional (`--positive_count`) | Always computed |
| Strand bias | Not computed | Fisher's exact test (read + fragment level) |
| Fractional depth | `--fragment_fractional_weight` | Not implemented |
| Parallelism | OpenMP block-based | Rayon per-variant |

## Related

- [Architecture](architecture.md) — System design overview
- [Input Formats](input-formats.md) — VCF and MAF specifications
- [Glossary](glossary.md) — Term definitions
- [CLI Run Command](../cli/run.md) — All parameter options
