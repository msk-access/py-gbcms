"""
Pipeline Orchestrator: Manages the execution flow of gbcms.

This module handles:
1. Reading variants from input (VCF/MAF).
2. Iterating over samples (BAM files).
3. Running the Rust-based counting engine for each sample.
4. Writing results to per-sample output files.
"""

from pathlib import Path
from typing import List

import gbcms_rs

from .models.core import GbcmsConfig, Variant, OutputFormat
from .io.input import VcfReader, MafReader, VariantReader
from .io.output import MafWriter, VcfWriter

class Pipeline:
    def __init__(self, config: GbcmsConfig):
        self.config = config
        
    def run(self):
        """Execute the pipeline."""
        print(f"Starting gbcms pipeline with config: {self.config}")
        
        # 1. Load Variants
        print(f"Loading variants from {self.config.variant_file}...")
        variants = self._load_variants()
        print(f"Loaded {len(variants)} variants.")
        
        if not variants:
            print("No variants found. Exiting.")
            return

        # 2. Prepare Rust Variants
        rs_variants = [
            gbcms_rs.Variant(v.chrom, v.pos, v.ref, v.alt, v.variant_type.value)
            for v in variants
        ]
        
        # 3. Process Each Sample
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        
        for sample_name, bam_path in self.config.bam_files.items():
            print(f"Processing sample: {sample_name} ({bam_path})...")
            
            try:
                # Run Rust Engine
                counts_list = gbcms_rs.count_bam(
                    str(bam_path),
                    rs_variants,
                    min_mapq=self.config.min_mapping_quality,
                    min_baseq=self.config.min_base_quality,
                    filter_duplicates=self.config.filter_duplicates,
                    filter_secondary=self.config.filter_secondary,
                    filter_supplementary=self.config.filter_supplementary
                )
                
                # Write Output
                self._write_output(sample_name, variants, counts_list)
                print(f"Finished sample: {sample_name}")
                
            except Exception as e:
                print(f"Error processing sample {sample_name}: {e}")
                # Continue to next sample? Or fail?
                # For now, log and continue.
                continue
                
        print("Pipeline completed.")

    def _load_variants(self) -> List[Variant]:
        """Load variants based on file extension."""
        path = self.config.variant_file
        reader: VariantReader
        
        if path.suffix.lower() in ['.vcf', '.gz']: # .vcf.gz handled by pysam
            reader = VcfReader(path)
        elif path.suffix.lower() == '.maf':
            reader = MafReader(path, fasta_path=self.config.fasta_file)
        else:
            # Default to VCF if unknown? Or error?
            # Let's assume VCF for now or check content.
            # But extension check is safer for now.
            raise ValueError(f"Unsupported variant file format: {path.suffix}")
            
        variants = list(reader)
        if hasattr(reader, 'close'):
            reader.close()
            
        return variants

    def _write_output(self, sample_name: str, variants: List[Variant], counts_list: List[gbcms_rs.BaseCounts]):
        """Write results to output file."""
        ext = "vcf" if self.config.output_format == OutputFormat.VCF else "maf"
        output_path = self.config.output_dir / f"{sample_name}.{ext}"
        
        writer = None
        if self.config.output_format == OutputFormat.VCF:
            writer = VcfWriter(output_path, sample_name=sample_name)
        else:
            writer = MafWriter(output_path)
            
        try:
            for variant, counts in zip(variants, counts_list):
                writer.write(variant, counts, sample_name=sample_name)
        finally:
            writer.close()
        
        print(f"Results written to {output_path}")
