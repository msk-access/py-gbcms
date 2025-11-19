"""
Output Writers: Formatting results for VCF and MAF.

This module provides classes to write processed variants and their counts
to output files, handling format-specific columns and headers.
"""

import csv
from pathlib import Path
from typing import List, Dict, Any

import pysam

from ..models.core import Variant, VariantType
# Import BaseCounts from Rust extension via a wrapper or direct import if possible
# For now, we assume we receive objects with the attributes defined in Rust
# We can define a Protocol for type hinting if needed.

class OutputWriter:
    """Abstract base class for output writers."""
    
    def write(self, variant: Variant, counts: Any):
        raise NotImplementedError
        
    def close(self):
        pass


class MafWriter(OutputWriter):
    """Writes results to a MAF-like file (Fillout format)."""
    
    def __init__(self, path: Path):
        self.path = path
        self.file = open(path, 'w')
        self.writer = None
        self._headers_written = False
        
    def _init_writer(self):
        # Standard GDC MAF columns (plus our custom ones)
        # Based on GDC MAF Format v1.0.0
        self.fieldnames = [
            "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", 
            "Start_Position", "End_Position", "Strand", "Variant_Classification", 
            "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", 
            "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode", 
            "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
            "Tumor_Validation_Allele1", "Tumor_Validation_Allele2", 
            "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2", 
            "Verification_Status", "Validation_Status", "Mutation_Status", 
            "Sequencing_Phase", "Sequence_Source", "Validation_Method", "Score", 
            "BAM_File", "Sequencer", "Tumor_Sample_UUID", "Matched_Norm_Sample_UUID",
            "HGVSc", "HGVSp", "HGVSp_Short", "Transcript_ID", "Exon_Number", 
            "t_depth", "t_ref_count", "t_alt_count", "n_depth", "n_ref_count", "n_alt_count", 
            "all_effects", "Allele", "Gene", "Feature", "Feature_type", "Consequence", 
            "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", 
            "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "SYMBOL", "SYMBOL_SOURCE", 
            "HGNC_ID", "BIOTYPE", "CANONICAL", "CCDS", "ENSP", "SWISSPROT", "TREMBL", "UNIPARC", 
            "RefSeq", "SIFT", "PolyPhen", "EXON", "INTRON", "DOMAINS", "GMAF", "AFR_MAF", 
            "AMR_MAF", "ASN_MAF", "EUR_MAF", "AA_MAF", "EA_MAF", "CLIN_SIG", "SOMATIC", 
            "PUBMED", "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", 
            "IMPACT", "PICK", "VARIANT_CLASS", "TSL", "HGVS_OFFSET", "PHENO", "MINIMISED", 
            "ExAC_AF", "ExAC_AF_AFR", "ExAC_AF_AMR", "ExAC_AF_EAS", "ExAC_AF_FIN", 
            "ExAC_AF_NFE", "ExAC_AF_OTH", "ExAC_AF_SAS", "GENE_PHENO", "FILTER", 
            "flanking_bps", "vcf_id", "vcf_qual", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", 
            "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF", 
            "gnomAD_SAS_AF", "vcf_pos", "vcf_region",
            # Custom columns
            "t_total_count", "t_vaf",
            "t_ref_count_fragment", "t_alt_count_fragment", "t_total_count_fragment", "t_vaf_fragment",
            "strand_bias_p_value", "strand_bias_odds_ratio",
            "fragment_strand_bias_p_value", "fragment_strand_bias_odds_ratio"
        ]
        self.writer = csv.DictWriter(self.file, fieldnames=self.fieldnames, delimiter='\t')
        self.writer.writeheader()
        self._headers_written = True
        
    def write(self, variant: Variant, counts: Any, sample_name: str = "TUMOR"):
        if not self._headers_written:
            self._init_writer()
            
        # Calculate VAFs
        total = counts.rd + counts.ad
        vaf = counts.ad / total if total > 0 else 0.0
        
        total_frag = counts.rdf + counts.adf
        vaf_frag = counts.adf / total_frag if total_frag > 0 else 0.0
        
        # MAF Coordinates (1-based)
        start_pos = variant.pos + 1
        end_pos = start_pos
        
        if variant.variant_type == "DELETION":
             end_pos = start_pos + len(variant.ref) - 1
        elif variant.variant_type == "INSERTION":
             # MAF for insertion: Start and End are the same (anchor), or Start=Anchor, End=Anchor+1?
             # GDC: Start_Position is the last base of the reference allele (anchor).
             # End_Position is Start_Position + 1.
             # Let's follow GDC convention if possible, or stick to VCF-like anchor.
             # For now, let's keep it simple: Start=End=Anchor for Ins?
             # Actually, standard MAF usually has Start=End for insertions (between bases).
             end_pos = start_pos + 1 # To indicate range?
        
        # Populate row with defaults for missing fields
        row = {k: "" for k in self.fieldnames}
        
        # Fill known fields
        row.update({
            "Chromosome": variant.chrom,
            "Start_Position": start_pos,
            "End_Position": end_pos,
            "Reference_Allele": variant.ref,
            "Tumor_Seq_Allele2": variant.alt,
            "Tumor_Sample_Barcode": sample_name,
            "Variant_Type": variant.variant_type,
            "t_ref_count": counts.rd,
            "t_alt_count": counts.ad,
            "t_total_count": counts.dp, 
            "t_vaf": f"{vaf:.4f}",
            "t_ref_count_fragment": counts.rdf,
            "t_alt_count_fragment": counts.adf,
            "t_total_count_fragment": counts.dpf,
            "t_vaf_fragment": f"{vaf_frag:.4f}",
            "strand_bias_p_value": f"{counts.sb_pval:.4e}",
            "strand_bias_odds_ratio": f"{counts.sb_or:.4f}",
            "fragment_strand_bias_p_value": f"{counts.fsb_pval:.4e}",
            "fragment_strand_bias_odds_ratio": f"{counts.fsb_or:.4f}",
            "vcf_region": f"{variant.chrom}:{start_pos}-{end_pos}", # Simple region string
            "vcf_pos": start_pos,
        })
        
        if variant.original_id:
            row["vcf_id"] = variant.original_id

        self.writer.writerow(row)
        
    def close(self):
        self.file.close()

class VcfWriter(OutputWriter):
    """Writes results to a VCF file."""
    
    def __init__(self, path: Path, sample_name: str = "SAMPLE"):
        self.path = path
        self.sample_name = sample_name
        self.file = open(path, 'w')
        self._headers_written = False
        
    def _write_header(self):
        # Minimal VCF header
        headers = [
            "##fileformat=VCFv4.2",
            "##source=gbcms_v2",
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
            '##INFO=<ID=SB_PVAL,Number=1,Type=Float,Description="Fisher strand bias p-value">',
            '##INFO=<ID=FSB_PVAL,Number=1,Type=Float,Description="Fisher fragment strand bias p-value">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">',
            '##FORMAT=<ID=RDF,Number=1,Type=Integer,Description="Ref Fragment Count">',
            '##FORMAT=<ID=ADF,Number=1,Type=Integer,Description="Alt Fragment Count">',
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.sample_name}"
        ]
        self.file.write("\n".join(headers) + "\n")
        self._headers_written = True
        
    def write(self, variant: Variant, counts: Any, sample_name: str = "SAMPLE"):
        if not self._headers_written:
            self._write_header()
            
        # VCF POS is 1-based
        pos = variant.pos + 1
        
        # INFO fields
        info = f"DP={counts.dp};SB_PVAL={counts.sb_pval:.4e};FSB_PVAL={counts.fsb_pval:.4e}"
        
        # FORMAT fields
        # GT: Simple 0/1 if alt > 0? Or ./1?
        # Let's assume 0/1 if we have alt counts, else 0/0
        gt = "0/1" if counts.ad > 0 else "0/0"
        ad = f"{counts.rd},{counts.ad}"
        dp = counts.dp
        rdf = counts.rdf
        adf = counts.adf
        
        format_str = "GT:AD:DP:RDF:ADF"
        sample_data = f"{gt}:{ad}:{dp}:{rdf}:{adf}"
        
        row = [
            variant.chrom,
            str(pos),
            variant.original_id or ".",
            variant.ref,
            variant.alt,
            ".", # QUAL
            ".", # FILTER
            info,
            format_str,
            sample_data
        ]
        
        self.file.write("\t".join(row) + "\n")
        
    def close(self):
        self.file.close()

# VCF Writer would require a template VCF or constructing a header from scratch.
# For MVP, MAF/Table output is often preferred for downstream analysis.
# We can add VCFWriter later if strictly needed, or use pysam.
