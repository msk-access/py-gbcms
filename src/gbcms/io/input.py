"""
Input Adapters: Handling VCF and MAF inputs.

This module provides classes to read variants from VCF and MAF files,
converting them into the internal normalized representation using CoordinateKernel.
"""

import csv
from pathlib import Path
from typing import Iterator, Optional

import pysam
from pydantic import ValidationError

from ..core.kernel import CoordinateKernel
from ..models.core import Variant, VariantType


class VariantReader:
    """Abstract base class for variant readers."""
    
    def __iter__(self) -> Iterator[Variant]:
        raise NotImplementedError


class VcfReader(VariantReader):
    """Reads variants from a VCF file."""
    
    def __init__(self, path: Path):
        self.path = path
        self._vcf = pysam.VariantFile(str(path))
        
    def __iter__(self) -> Iterator[Variant]:
        for record in self._vcf:
            # VCF coordinates are 1-based
            # pysam converts them to 0-based automatically?
            # pysam.VariantFile returns 0-based pos (start)
            # BUT CoordinateKernel.vcf_to_internal expects 1-based VCF POS.
            # Let's check pysam documentation or behavior.
            # pysam record.pos is 0-based. record.start is 0-based.
            # The VCF file itself has 1-based POS.
            # If we use record.pos + 1, we get the VCF POS.
            
            # Handle multiple ALTs
            for alt in record.alts or []:
                # VCF POS is record.pos + 1
                yield CoordinateKernel.vcf_to_internal(
                    chrom=record.chrom,
                    pos=record.pos + 1, 
                    ref=record.ref,
                    alt=alt,
                    original_id=record.id
                )
                
    def close(self):
        self._vcf.close()


class MafReader(VariantReader):
    """Reads variants from a MAF file."""
    
    def __init__(self, path: Path, fasta_path: Optional[Path] = None):
        self.path = path
        self.fasta = pysam.FastaFile(str(fasta_path)) if fasta_path else None
        
    def __iter__(self) -> Iterator[Variant]:
        with open(self.path, 'r') as f:
            # Skip comments
            while True:
                pos = f.tell()
                line = f.readline()
                if not line.startswith('#'):
                    f.seek(pos)
                    break
            
            reader = csv.DictReader(f, delimiter='\t')
            
            for row in reader:
                try:
                    chrom = row['Chromosome']
                    start_pos = int(row['Start_Position'])
                    ref = row['Reference_Allele']
                    alt = row['Tumor_Seq_Allele2'] # Standard MAF alt column
                    
                    # Normalize Indels if FASTA is available
                    if self.fasta and (ref == '-' or alt == '-'):
                        if ref == '-': # Insertion
                            # MAF Start_Position is the base BEFORE the insertion (anchor)
                            # 1-based coordinate
                            anchor_pos_1based = start_pos
                            anchor_pos_0based = anchor_pos_1based - 1
                            
                            # Fetch anchor base
                            # Try normalized and original chromosome names
                            norm_chrom = CoordinateKernel.normalize_chromosome(chrom)
                            try:
                                anchor_base = self.fasta.fetch(norm_chrom, anchor_pos_0based, anchor_pos_0based + 1).upper()
                            except (KeyError, ValueError):
                                try:
                                    anchor_base = self.fasta.fetch(chrom, anchor_pos_0based, anchor_pos_0based + 1).upper()
                                except (KeyError, ValueError):
                                    # If both fail, we can't normalize. Skip or raise?
                                    # For now, skip/log
                                    continue
                            
                            # VCF Style:
                            # POS = anchor_pos_1based
                            # REF = anchor_base
                            # ALT = anchor_base + inserted_seq
                            vcf_pos = anchor_pos_1based
                            vcf_ref = anchor_base
                            vcf_alt = anchor_base + alt
                            
                        else: # Deletion (alt == '-')
                            # MAF Start_Position is the FIRST DELETED base
                            # Anchor is the base before that
                            first_deleted_1based = start_pos
                            anchor_pos_1based = first_deleted_1based - 1
                            anchor_pos_0based = anchor_pos_1based - 1
                            
                            # Fetch anchor base
                            norm_chrom = CoordinateKernel.normalize_chromosome(chrom)
                            try:
                                anchor_base = self.fasta.fetch(norm_chrom, anchor_pos_0based, anchor_pos_0based + 1).upper()
                            except (KeyError, ValueError):
                                try:
                                    anchor_base = self.fasta.fetch(chrom, anchor_pos_0based, anchor_pos_0based + 1).upper()
                                except (KeyError, ValueError):
                                    continue
                            
                            # VCF Style:
                            # POS = anchor_pos_1based
                            # REF = anchor_base + deleted_seq
                            # ALT = anchor_base
                            vcf_pos = anchor_pos_1based
                            vcf_ref = anchor_base + ref
                            vcf_alt = anchor_base
                            
                        yield CoordinateKernel.vcf_to_internal(
                            chrom=chrom,
                            pos=vcf_pos,
                            ref=vcf_ref,
                            alt=vcf_alt
                        )
                    else:
                        # Fallback to old behavior or direct mapping for SNPs
                        # For SNPs, MAF Start_Position == VCF POS
                        if len(ref) == len(alt) == 1 and ref != '-' and alt != '-':
                             yield CoordinateKernel.vcf_to_internal(
                                chrom=chrom,
                                pos=start_pos,
                                ref=ref,
                                alt=alt
                            )
                        else:
                             # Fallback for complex/unhandled without FASTA
                             # This might fail in Rust engine if it expects anchor
                             yield CoordinateKernel.maf_to_internal(
                                chrom=chrom,
                                start_pos=start_pos,
                                end_pos=int(row['End_Position']),
                                ref=ref,
                                alt=alt
                            )

                except (KeyError, ValueError, ValidationError) as e:
                    # Log warning or skip malformed lines
                    continue
                    
    def close(self):
        if self.fasta:
            self.fasta.close()

class ReferenceChecker:
    """
    Utility to check variants against a reference FASTA.
    Ensures that the REF allele matches the genome.
    """
    
    def __init__(self, fasta_path: Path):
        self.fasta = pysam.FastaFile(str(fasta_path))
        
    def validate(self, variant: Variant) -> bool:
        """
        Check if variant REF matches reference genome.
        """
        # Variant pos is 0-based.
        # Fetch sequence of length REF
        try:
            ref_seq = self.fasta.fetch(
                variant.chrom, 
                variant.pos, 
                variant.pos + len(variant.ref)
            )
            return ref_seq.upper() == variant.ref.upper()
        except (ValueError, KeyError):
            return False
            
    def close(self):
        self.fasta.close()
