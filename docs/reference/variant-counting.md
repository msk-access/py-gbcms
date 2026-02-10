# Variant Counting

How py-gbcms counts reads for each variant type — with visual examples.

## Overview

For every variant, py-gbcms fetches reads overlapping the variant position from each BAM file, applies [read filters](#read-filters), and then dispatches to a **type-specific** allele checker. Each read is classified as supporting the **reference** allele, the **alternate** allele, or **neither**.

```mermaid
flowchart LR
    Fetch([📥 Fetch Reads]):::fetch --> Filter([🔍 Apply Filters]):::filter
    Filter --> Dispatch{Variant Type?}
    Dispatch -->|SNP| SNP[check_snp]
    Dispatch -->|Insertion| INS[check_insertion]
    Dispatch -->|Deletion| DEL[check_deletion]
    Dispatch -->|MNP| MNP[check_mnp]
    Dispatch -->|Complex| CPX[check_complex]
    SNP --> Count([📊 Update Counts]):::count
    INS --> Count
    DEL --> Count
    MNP --> Count
    CPX --> Count

    classDef fetch fill:#4a90d9,color:#fff,stroke:#2c6fad,stroke-width:2px;
    classDef filter fill:#e67e22,color:#fff,stroke:#bf6516,stroke-width:2px;
    classDef count fill:#27ae60,color:#fff,stroke:#1e8449,stroke-width:2px;
```

---

## Read Filter Pipeline

Before any allele checking begins, every read passes through a **filter cascade**. Reads that fail any enabled filter are discarded. The order matches the Rust engine implementation.

```mermaid
flowchart TD
    Start([📖 Read from BAM]):::process

    subgraph Filters [Quality & Alignment Filters]
        direction TB
        F1{Duplicate?}:::decision
        F2{Secondary?}:::decision
        F3{Supplementary?}:::decision
        F4{QC Failed?}:::decision
        F5{Improper pair?}:::decision
        F6{Contains indel?}:::decision
        F7{MAPQ ≥ 20?}:::decision
    end

    Drop{{❌ Discard Read}}:::discard
    Pass([✅ Pass to Allele Checker]):::success

    Start --> F1

    F1 -- "Yes" --> Drop
    F1 -- "No" --> F2

    F2 -- "Yes" --> Drop
    F2 -- "No" --> F3

    F3 -- "Yes" --> Drop
    F3 -- "No" --> F4

    F4 -- "Yes" --> Drop
    F4 -- "No" --> F5

    F5 -- "Yes" --> Drop
    F5 -- "No" --> F6

    F6 -- "Yes" --> Drop
    F6 -- "No" --> F7

    F7 -- "Below 20" --> Drop
    F7 -- "Pass" --> Pass

    classDef process fill:#9b59b6,color:#fff,stroke:#7d3c98,stroke-width:2px;
    classDef decision fill:#f39c12,color:#fff,stroke:#d68910,stroke-width:2px;
    classDef discard fill:#e74c3c,color:#fff,stroke:#c0392b,stroke-width:2px;
    classDef success fill:#27ae60,color:#fff,stroke:#1e8449,stroke-width:2px;
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
    Start([🧬 SNP Check]):::start --> Walk[Walk CIGAR to variant pos]
    Walk --> Found{Position found?}
    Found -->|No| Neither1([Neither]):::neither
    Found -->|Yes| BQ{Base quality ≥ min_baseq?}
    BQ -->|No| Neither2([Neither]):::neither
    BQ -->|Yes| Compare[Compare base to REF and ALT]
    Compare --> IsRef{base == REF?}
    IsRef -->|Yes| Ref([✅ REF]):::ref
    IsRef -->|No| IsAlt{base == ALT?}
    IsAlt -->|Yes| Alt([🔴 ALT]):::alt
    IsAlt -->|No| Neither3([Neither]):::neither

    classDef start fill:#9b59b6,color:#fff,stroke:#7d3c98,stroke-width:2px;
    classDef ref fill:#27ae60,color:#fff,stroke:#1e8449,stroke-width:2px;
    classDef alt fill:#e74c3c,color:#fff,stroke:#c0392b,stroke-width:2px;
    classDef neither fill:#95a5a6,color:#fff,stroke:#7f8c8d,stroke-width:2px;
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

The insertion check uses a **single CIGAR walk** with two strategies: strict match (fast path) and windowed scan (±5bp fallback).

```mermaid
flowchart TD
    Start([🧬 Insertion Check]):::start --> Walk[Walk CIGAR left → right]
    Walk --> Found{Match block contains anchor?}
    Found -->|No| WindowCheck
    Found -->|Yes| AtEnd{Anchor at end of block?}
    AtEnd -->|No| RefCov[Mark ref coverage]
    AtEnd -->|Yes| NextOp{Next op is Ins?}
    NextOp -->|No| RefCov
    NextOp -->|Yes| Strict{Length + seq match?}
    Strict -->|Yes| StrictAlt([🔴 ALT - strict]):::alt
    Strict -->|No| RefCov
    RefCov --> WindowCheck

    WindowCheck{Ins within ±5bp window?}
    WindowCheck -->|No| Continue[Continue CIGAR walk]
    WindowCheck -->|Yes| S1{S1: Seq matches?}
    S1 -->|No| Continue
    S1 -->|Yes| S3{S3: Anchor base matches ref?}
    S3 -->|No| Continue
    S3 -->|Yes| S2[S2: Track closest match]
    S2 --> Continue

    Continue --> MoreOps{More CIGAR ops?}
    MoreOps -->|Yes| Walk
    MoreOps -->|No| Eval{Windowed candidate found?}
    Eval -->|Yes| WinAlt([🔴 ALT - windowed]):::alt
    Eval -->|No| HasRef{Read covered anchor?}
    HasRef -->|Yes| Ref([✅ REF]):::ref
    HasRef -->|No| Neither([Neither]):::neither

    classDef start fill:#9b59b6,color:#fff,stroke:#7d3c98,stroke-width:2px;
    classDef ref fill:#27ae60,color:#fff,stroke:#1e8449,stroke-width:2px;
    classDef alt fill:#e74c3c,color:#fff,stroke:#c0392b,stroke-width:2px;
    classDef neither fill:#95a5a6,color:#fff,stroke:#7f8c8d,stroke-width:2px;
```

#### Visual Example

```
Variant: chr1:100 A→ATG (insertion of TG after anchor A)

Reference:     5'─ ...C  G  A ── C  G  T  A... ─3'
                          99 100  101 102
                               ▲
                          anchor pos

Read 1 (ALT, strict):  CIGAR = 5M 2I 5M
               5'─ ...C  G  A [T  G] C  G  T  A... ─3'
                               └──┘
                          inserted bases at anchor → ALT ✅ (strict)

Read 2 (ALT, windowed): CIGAR = 7M 2I 3M
               5'─ ...C  G  A  C  G [T  G] T  A... ─3'
                                          └──┘
                    insertion shifted +2bp, same seq → ALT ✅ (windowed)

Read 3 (REF):  CIGAR = 10M
               5'─ ...C  G  A  C  G  T  A... ─3'
                               ↑
                          no insertion after anchor → REF ✅

Read 4 (wrong seq): CIGAR = 5M 2I 5M
               5'─ ...C  G  A [C  C] C  G  T  A... ─3'
                               └──┘
                          insertion bases ≠ expected TG → S1 reject
```

!!! tip "Why anchor-based?"
    In VCF format, insertions are represented as `REF=A, ALT=ATG` where the first base `A` is the anchor — it's not part of the inserted sequence. py-gbcms uses this anchor to locate where the insertion should occur in the CIGAR string.

!!! info "MAF → VCF Anchor Conversion"
    MAF format represents indels differently — using `-` dashes instead of an anchor base. py-gbcms automatically converts MAF indels to VCF-style anchor-based representation at input time, which requires a **reference FASTA** (`--fasta`).

    **Insertion** (e.g., insert `TG` after chr1:100):

    | | MAF | VCF (internal) |
    |:--|:--|:--|
    | Position | `Start_Position=100` (anchor) | `POS=100` |
    | REF | `-` | `A` (fetched from FASTA) |
    | ALT | `TG` | `ATG` (anchor + inserted seq) |

    **Deletion** (e.g., delete `CG` at chr1:101–102):

    | | MAF | VCF (internal) |
    |:--|:--|:--|
    | Position | `Start_Position=101` (first deleted base) | `POS=100` (anchor = base before) |
    | REF | `CG` | `ACG` (anchor + deleted seq) |
    | ALT | `-` | `A` (anchor only) |

    Note that for **deletions**, MAF `Start_Position` points to the *first deleted base*, not the anchor. py-gbcms shifts back by one position and fetches the preceding base from the FASTA as the anchor.

    See [Input Formats](input-formats.md#maf-indel-normalization) for full details.

---

### Deletion

Bases deleted after an **anchor** position. The logic mirrors insertion but looks for `Del` CIGAR operations.

| Property | Value |
|:---------|:------|
| Detection | `len(REF) > 1 && len(ALT) == 1` |
| Position | 0-based index of the **anchor** base |
| Quality check | None (CIGAR-structural check only) |

#### Algorithm Flow

Same single-walk strategy as insertion, with Safeguard 3 checking deleted reference bases.

```mermaid
flowchart TD
    Start([🧬 Deletion Check]):::start --> Walk[Walk CIGAR left → right]
    Walk --> Found{Match block contains anchor?}
    Found -->|No| WindowCheck
    Found -->|Yes| AtEnd{Anchor at end of block?}
    AtEnd -->|No| RefCov[Mark ref coverage]
    AtEnd -->|Yes| NextOp{Next op is Del?}
    NextOp -->|No| RefCov
    NextOp -->|Yes| Strict{Length matches?}
    Strict -->|Yes| StrictAlt([🔴 ALT - strict]):::alt
    Strict -->|No| RefCov
    RefCov --> WindowCheck

    WindowCheck{Del within ±5bp window?}
    WindowCheck -->|No| Continue[Continue CIGAR walk]
    WindowCheck -->|Yes| S1{S1: Del length matches?}
    S1 -->|No| Continue
    S1 -->|Yes| S3{S3: Ref bases at shift match?}
    S3 -->|No| Continue
    S3 -->|Yes| S2[S2: Track closest match]
    S2 --> Continue

    Continue --> MoreOps{More CIGAR ops?}
    MoreOps -->|Yes| Walk
    MoreOps -->|No| Eval{Windowed candidate found?}
    Eval -->|Yes| WinAlt([🔴 ALT - windowed]):::alt
    Eval -->|No| HasRef{Read covered anchor?}
    HasRef -->|Yes| Ref([✅ REF]):::ref
    HasRef -->|No| Neither([Neither]):::neither

    classDef start fill:#9b59b6,color:#fff,stroke:#7d3c98,stroke-width:2px;
    classDef ref fill:#27ae60,color:#fff,stroke:#1e8449,stroke-width:2px;
    classDef alt fill:#e74c3c,color:#fff,stroke:#c0392b,stroke-width:2px;
    classDef neither fill:#95a5a6,color:#fff,stroke:#7f8c8d,stroke-width:2px;
```

#### Visual Example

```
Variant: chr1:100 ACG→A (deletion of CG after anchor A)

Reference:     5'─ ...T  G  A  C  G  T  A  C... ─3'
                          99 100 101 102
                               ▲
                          anchor pos

Read 1 (ALT, strict):  CIGAR = 5M 2D 5M
               5'─ ...T  G  A  ──  ──  T  A  C... ─3'
                               └─────┘
                          2bp deletion at anchor → ALT ✅ (strict)

Read 2 (ALT, windowed): CIGAR = 7M 2D 3M
               5'─ ...T  G  A  C  G  T  ──  ──  A  C... ─3'
                                              └─────┘
                  deletion shifted +2bp, same length, ref matches → ALT ✅ (windowed)
                  (S3 verifies ref bases at shifted pos match deleted 'CG')

Read 3 (REF):  CIGAR = 12M
               5'─ ...T  G  A  C  G  T  A  C... ─3'
                               ↑
                          no deletion after anchor → REF ✅

Read 4 (wrong length): CIGAR = 5M 3D 5M
               5'─ ...T  G  A  ──  ──  ──  A  C... ─3'
                               └────────┘
                          3bp deletion ≠ expected 2bp → S1 reject
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
    Start([🧬 MNP Check]):::start --> Find[Find read position of first base]
    Find --> Found{Position found?}
    Found -->|No| Neither1([Neither]):::neither
    Found -->|Yes| Cover{Read covers entire MNP region?}
    Cover -->|No| Neither2([Neither]):::neither
    Cover -->|Yes| Loop[For each position in MNP region]

    subgraph PerBase [Per-Base Validation]
        direction TB
        BQ{Base quality ≥ min_baseq?}
        Track[Compare base to REF and ALT]
        More{More positions?}
    end

    Loop --> BQ
    BQ -->|No| Neither3([Neither]):::neither
    BQ -->|Yes| Track
    Track --> More
    More -->|Yes| BQ
    More -->|No| Contig{No indels within MNP region?}

    Contig -->|Indel found| Neither4([Neither]):::neither
    Contig -->|Contiguous| Final{All bases match?}
    Final -->|All match ALT| Alt([🔴 ALT]):::alt
    Final -->|All match REF| Ref([✅ REF]):::ref
    Final -->|Mixed| Neither5([Neither]):::neither

    classDef start fill:#9b59b6,color:#fff,stroke:#7d3c98,stroke-width:2px;
    classDef ref fill:#27ae60,color:#fff,stroke:#1e8449,stroke-width:2px;
    classDef alt fill:#e74c3c,color:#fff,stroke:#c0392b,stroke-width:2px;
    classDef neither fill:#95a5a6,color:#fff,stroke:#7f8c8d,stroke-width:2px;
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

Variants where REF and ALT differ in both sequence **and** length. This is the catch-all category that uses a sophisticated **haplotype reconstruction** algorithm with **quality-aware masked comparison**.

| Property | Value |
|:---------|:------|
| Detection | Fallback for all other combinations |
| Position | 0-based index of the first reference base |
| Quality check | **Masked comparison** — bases below `--min-baseq` are masked out and cannot vote for either allele |

#### Algorithm Flow — Haplotype Reconstruction + Masked Comparison

```mermaid
flowchart TD
    Start([🧬 Complex Check]):::start --> Region["Define region: pos .. pos+len REF"]:::info
    Region --> Init["Init reconstructed_seq + quals_per_base"]

    subgraph CIGARWalk [Phase 1: CIGAR Walk — Haplotype Reconstruction]
        direction TB
        Walk[Walk each CIGAR op]
        OpType{CIGAR op type?}
        Walk --> OpType
        OpType -->|"M / = / X"| Match["Append bases + quals"]
        OpType -->|"I"| InsCheck["Append inserted bases + quals"]
        OpType -->|"D / N"| AdvRef[Advance ref_pos only]
        OpType -->|"S"| AdvRead[Advance read_pos only]
        OpType -->|"H / P"| Skip[No action]
        Match --> Next[Next op]
        InsCheck --> Next
        AdvRef --> Next
        AdvRead --> Next
        Skip --> Next
        Next --> More{More CIGAR ops?}
        More -->|Yes| Walk
    end

    Init --> Walk
    More -->|No| LenCheck

    subgraph MaskedCompare ["Phase 2: Masked Comparison (Reliable Intersection)"]
        direction TB
        LenCheck{Check recon length}
        LenCheck -->|"== ALT == REF"| CaseA["Case A: Dual masked compare"]
        LenCheck -->|"== ALT only"| CaseB["Case B: ALT-only masked compare"]
        LenCheck -->|"== REF only"| CaseC["Case C: REF-only masked compare"]
        LenCheck -->|Neither| Neither4([Neither]):::neither

        CaseA --> Mask1["Mask bases with qual < min_baseq"]
        Mask1 --> Reliable1{reliable_count > 0?}
        Reliable1 -->|No| Neither5([Neither]):::neither
        Reliable1 -->|Yes| Ambig{Matches BOTH on reliable?}
        Ambig -->|Yes| Neither6(["Neither (ambiguous)"]):::neither
        Ambig -->|No| AltOnly{Matches ALT only?}
        AltOnly -->|Yes| Alt1([🔴 ALT]):::alt
        AltOnly -->|No| RefOnly{Matches REF only?}
        RefOnly -->|Yes| Ref1([✅ REF]):::ref
        RefOnly -->|No| Neither7([Neither]):::neither

        CaseB --> Mask2["Mask low-qual bases"]
        Mask2 --> ReliableB{reliable > 0 AND 0 mismatches?}
        ReliableB -->|Yes| Alt2([🔴 ALT]):::alt
        ReliableB -->|No| Neither8([Neither]):::neither

        CaseC --> Mask3["Mask low-qual bases"]
        Mask3 --> ReliableC{reliable > 0 AND 0 mismatches?}
        ReliableC -->|Yes| Ref2([✅ REF]):::ref
        ReliableC -->|No| Neither9([Neither]):::neither
    end

    classDef start fill:#9b59b6,color:#fff,stroke:#7d3c98,stroke-width:2px;
    classDef info fill:#3498db,color:#fff,stroke:#2471a3,stroke-width:2px;
    classDef ref fill:#27ae60,color:#fff,stroke:#1e8449,stroke-width:2px;
    classDef alt fill:#e74c3c,color:#fff,stroke:#c0392b,stroke-width:2px;
    classDef neither fill:#95a5a6,color:#fff,stroke:#7f8c8d,stroke-width:2px;
```

!!! info "Masked Comparison — Why Not Exact Match?"
    A single sequencing error at a low-quality base would discard the entire read
    under exact matching. The masked approach ignores unreliable bases, allowing the
    remaining high-quality bases to determine allele support. When masking causes both
    alleles to match (ambiguity), the read is safely discarded rather than guessed.

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
    subgraph ReadLevel [📖 Read-Level Counting]
        R1[Read 1 fwd → ALT]:::altread
        R2[Read 2 rev → ALT]:::altread
        R3[Read 3 fwd → REF]:::refread
        R4[Read 4 rev → REF]:::refread
        R5[Read 5 fwd → ALT]:::altread
        R6[Read 6 rev → REF]:::refread
    end

    subgraph FragLevel [🧬 Fragment-Level Counting]
        direction TB
        F1[Fragment A: R1+R2 → ALT]:::altfrag
        F2[Fragment B: R3+R4 → REF]:::reffrag
        F3[Fragment C: R5 only → ALT]:::altfrag
        F4[Fragment D: R6 only → REF]:::reffrag
    end

    R1 --> F1
    R2 --> F1
    R3 --> F2
    R4 --> F2
    R5 --> F3
    R6 --> F4

    subgraph Results [📊 Final Counts]
        ReadCounts[Reads: DP=6  RD=3  AD=3]
        FragCounts[Fragments: DPF=4  RDF=2  ADF=2]
    end

    F1 --> FragCounts
    F2 --> FragCounts
    F3 --> FragCounts
    F4 --> FragCounts

    classDef altread fill:#e74c3c15,stroke:#e74c3c;
    classDef refread fill:#27ae6015,stroke:#27ae60;
    classDef altfrag fill:#9b59b615,stroke:#9b59b6;
    classDef reffrag fill:#27ae6015,stroke:#27ae60;
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

Fragment counting collapses read pairs into a single observation per fragment using **quality-weighted consensus**. Each fragment is tracked internally as a `FragmentEvidence` struct that records the best base quality seen for REF and ALT across both reads.

```mermaid
flowchart TD
    Start([For each read passing filters]):::start
    Start --> Hash["Hash QNAME → u64 key"]
    Hash --> Observe["Record allele + base quality into FragmentEvidence"]
    Observe --> Done[After all reads processed]
    Done --> ForEach[For each unique fragment]
    ForEach --> HasBoth{Both REF and ALT evidence?}
    HasBoth -->|No| Direct[Assign to whichever allele was seen]
    HasBoth -->|Yes| QualCheck{"Quality difference > threshold?"}
    QualCheck -->|"REF qual >> ALT qual"| Ref([✅ Count as REF]):::ref
    QualCheck -->|"ALT qual >> REF qual"| Alt([🔴 Count as ALT]):::alt
    QualCheck -->|"Within threshold"| Discard([⚪ Discard — ambiguous]):::discard
    Direct --> Count
    Ref --> Count
    Alt --> Count
    Discard --> DPF["Still counted in DPF"]:::info
    Count[Update RDF / ADF / DPF + strand counts]

    classDef start fill:#3498db,color:#fff,stroke:#2471a3,stroke-width:2px;
    classDef ref fill:#27ae60,color:#fff,stroke:#1e8449,stroke-width:2px;
    classDef alt fill:#e74c3c,color:#fff,stroke:#c0392b,stroke-width:2px;
    classDef discard fill:#95a5a6,color:#fff,stroke:#7f8c8d,stroke-width:2px;
    classDef info fill:#3498db15,stroke:#3498db;
```

#### Quality-Weighted Consensus

When R1 and R2 of a fragment **disagree** (one supports REF, the other ALT), the engine resolves the conflict using base quality scores:

| Scenario | Condition | Result |
|:---------|:----------|:-------|
| **REF wins** | `best_ref_qual > best_alt_qual + threshold` | Count as REF |
| **ALT wins** | `best_alt_qual > best_ref_qual + threshold` | Count as ALT |
| **Ambiguous** | Quality difference ≤ threshold | **Discard** (neither REF nor ALT) |

The threshold is configurable via `--fragment-qual-threshold` (default: **10**).

!!! important "Why discard instead of defaulting to REF?"
    Assigning ambiguous fragments to REF would systematically **deflate VAF** by inflating the reference count. In cfDNA sequencing where true variants can be at 0.1–1% VAF, this bias could mask real mutations. Discarding preserves an unbiased VAF estimate at the cost of slightly reduced power (fewer counted fragments).

!!! tip "Quality metric: DPF − (RDF + ADF)"
    Discarded fragments are still counted in **DPF** (total fragment depth) but *not* in RDF or ADF. The gap `DPF − (RDF + ADF)` reveals how many ambiguous fragments exist at a locus — a useful quality signal. A high gap suggests a noisy or error-prone site.

| Metric | Description |
|:-------|:------------|
| **DPF** | Fragment depth (all fragments, including discarded) |
| **RDF** / **ADF** | Reference / Alternate fragment counts (resolved only) |
| **RDF_fwd** / **RDF_rev** | Strand-specific reference fragment counts |
| **ADF_fwd** / **ADF_rev** | Strand-specific alternate fragment counts |

!!! tip "Why fragment counting matters for cfDNA"
    In cfDNA sequencing (like MSK-ACCESS), the same DNA fragment is often sequenced from both ends. Counting reads would double-count each fragment. Fragment-level counting gives a more accurate picture of the number of **unique molecules** supporting each allele.

---

### Strand Bias (Fisher's Exact Test)

Strand bias detects when an allele is disproportionately supported by reads on one strand — a common sequencing artifact.

```mermaid
flowchart LR
    subgraph Table [2×2 Contingency Table]
        direction TB
        Header["| · | Fwd | Rev |"]
        RefRow["| REF | a | b |"]
        AltRow["| ALT | c | d |"]
    end

    Table --> Fisher([Fisher Exact Test]):::test
    Fisher --> Pval[SB_pval]
    Fisher --> OR[SB_OR]

    Pval --> Interpret{p < 0.05?}
    Interpret -->|Yes| Bias{{⚠️ Possible Strand Bias}}:::warn
    Interpret -->|No| NoBias([✅ No Strand Bias]):::pass

    classDef test fill:#3498db,color:#fff,stroke:#2471a3,stroke-width:2px;
    classDef warn fill:#e74c3c,color:#fff,stroke:#c0392b,stroke-width:2px;
    classDef pass fill:#27ae60,color:#fff,stroke:#1e8449,stroke-width:2px;
```

Computed at **both** levels:

| Metric | Level | Description |
|:-------|:------|:------------|
| **SB_pval** / **SB_OR** | Read | Strand bias from individual reads |
| **FSB_pval** / **FSB_OR** | Fragment | Strand bias from collapsed fragments |

!!! example "Strand bias example"
    If a variant has `AD_fwd=15, AD_rev=1`, that's suspicious — almost all ALT-supporting reads are on the forward strand. Fisher's test would yield a low p-value, flagging this as a potential artifact rather than a true variant.

---

## Related

- [Architecture](architecture.md) — System design, full pipeline diagram, and comparison with original GBCMS
- [Input Formats](input-formats.md) — VCF and MAF specifications
- [Glossary](glossary.md) — Term definitions
- [CLI Run Command](../cli/run.md) — All parameter options
