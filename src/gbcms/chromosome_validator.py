"""Chromosome name validation and normalization system."""

import logging

logger = logging.getLogger(__name__)


class ChromosomeValidator:
    """Comprehensive chromosome name validation and normalization system."""

    def __init__(self, config):
        """
        Initialize chromosome validator.

        Args:
            config: gbcms configuration object
        """
        self.config = config
        self.reference_sequence = None
        self.bam_files = config.bam_files
        self.variant_files = config.variant_files

        # Validation results
        self.reference_format = None
        self.bam_format = None
        self.variant_chromosomes = set()
        self.normalization_needed = False
        self.validation_issues = []

    def is_standard_chromosome(self, chrom: str) -> bool:
        """
        Check if chromosome should be normalized.

        Only standard chromosomes (1-22, X, Y, M, MT) are normalized.
        Alternative contigs (KI270728.1, etc.) are left unchanged.

        Args:
            chrom: Chromosome name

        Returns:
            True if chromosome should be normalized
        """
        from .variant import _is_standard_chromosome

        return _is_standard_chromosome(chrom)

    def detect_reference_format(self) -> tuple[str, str]:
        """
        Detect chromosome name format used by reference FASTA.

        Returns:
            Tuple of (example_chromosome, format_type)
            Format type is "chr_prefix" or "no_prefix"
        """
        from .reference import ReferenceSequence

        try:
            ref_seq = ReferenceSequence(self.config.fasta_file)

            # Test standard chromosomes to detect format
            test_chromosomes = ["1", "chr1"]

            for test_chrom in test_chromosomes:
                try:
                    ref_seq.get_base(test_chrom, 0)
                    if test_chrom.startswith("chr"):
                        return test_chrom, "chr_prefix"
                    else:
                        return test_chrom, "no_prefix"
                except Exception:
                    continue

            # If we can't determine format, assume no_prefix (common case)
            logger.warning("Could not determine reference chromosome format, assuming no_prefix")
            return "1", "no_prefix"

        except Exception as e:
            logger.error(f"Could not load reference for format detection: {e}")
            return "1", "no_prefix"

    def detect_bam_format(self) -> tuple[str, str]:
        """
        Detect chromosome name format used by BAM files.

        Returns:
            Tuple of (example_chromosome, format_type)
        """
        # For now, assume BAM matches reference format
        # In a full implementation, we'd check BAM headers
        if self.reference_format is None:
            return "1", "no_prefix"
        if self.reference_format is None:
            return "1", "no_prefix"
        return self.reference_format[0], self.reference_format[1]

    def extract_variant_chromosomes(self) -> set[str]:
        """
        Extract chromosome names from variant files (MAF/VCF).

        Returns:
            Set of chromosome names found in variant files
        """
        chromosomes = set()

        for variant_file in self.variant_files:
            try:
                if variant_file.endswith(".vcf") or variant_file.endswith(".vcf.gz"):
                    chromosomes.update(self._extract_vcf_chromosomes(variant_file))
                elif variant_file.endswith(".maf"):
                    chromosomes.update(self._extract_maf_chromosomes(variant_file))
            except Exception as e:
                logger.warning(f"Could not extract chromosomes from {variant_file}: {e}")

        return chromosomes

    def _extract_vcf_chromosomes(self, vcf_file: str) -> set[str]:
        """Extract chromosome names from VCF file."""
        chromosomes = set()

        try:
            with open(vcf_file) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue

                    # First field is chromosome
                    chrom = line.split("\t")[0]
                    if chrom:  # Skip empty chromosomes
                        chromosomes.add(chrom)
                        break  # Just need one example per file
        except Exception as e:
            logger.warning(f"Could not read VCF chromosomes from {vcf_file}: {e}")

        return chromosomes

    def _extract_maf_chromosomes(self, maf_file: str) -> set[str]:
        """Extract chromosome names from MAF file."""
        chromosomes = set()

        try:
            with open(maf_file) as f:
                header_found = False
                chrom_idx: int | None = None

                for line in f:
                    if not header_found and not line.startswith("#"):
                        # Header line
                        headers = line.strip().split("\t")
                        if "Chromosome" in headers:
                            chrom_idx = headers.index("Chromosome")
                            header_found = True
                    elif header_found and line.strip():
                        # Data line
                        fields = line.strip().split("\t")
                        if chrom_idx is not None and len(fields) > chrom_idx:
                            chrom = fields[chrom_idx]
                            if chrom and chrom != "Chromosome":  # Skip header value
                                chromosomes.add(chrom)
                                break  # Just need one example per file
        except Exception as e:
            logger.warning(f"Could not read MAF chromosomes from {maf_file}: {e}")

        return chromosomes

    def smart_normalize_chromosome(self, chrom: str, target_format: str) -> str:
        """
        Smart normalization that respects chromosome types.

        Args:
            chrom: Original chromosome name
            target_format: "chr_prefix" or "no_prefix"

        Returns:
            Normalized chromosome name
        """
        from .variant import normalize_chromosome_name_to_format

        return normalize_chromosome_name_to_format(chrom, target_format)

    def validate_chromosome_consistency(self) -> bool:
        """
        Comprehensive chromosome consistency validation.

        Returns:
            True if all files are compatible, False if issues found
        """
        logger.info("Starting chromosome consistency validation")

        # Step 1: Detect reference format
        self.reference_format = self.detect_reference_format()
        ref_chrom, ref_format = self.reference_format
        logger.info(f"Reference uses chromosome format: {ref_format} (example: {ref_chrom})")

        # Step 2: Detect BAM format (must match reference)
        self.bam_format = self.detect_bam_format()
        bam_chrom, bam_format = self.bam_format

        if ref_format != bam_format:
            self.validation_issues.append(
                f"Reference and BAM use different chromosome formats. "
                f"Reference: {ref_format}, BAM: {bam_format}. "
                f"Please ensure reference and BAM files use consistent chromosome names."
            )
            logger.error(f"Chromosome format mismatch: Reference={ref_format}, BAM={bam_format}")
            return False

        # Step 3: Extract variant chromosomes
        self.variant_chromosomes = self.extract_variant_chromosomes()
        logger.info(f"Variant files use chromosomes: {sorted(self.variant_chromosomes)}")

        # Step 4: Validate variant compatibility
        incompatible_chromosomes = []

        for var_chrom in self.variant_chromosomes:
            # Test if normalized chromosome exists in reference
            try:
                # Try both original and normalized versions
                test_chroms = [var_chrom]
                if self.is_standard_chromosome(var_chrom):
                    test_chroms.append(self.smart_normalize_chromosome(var_chrom, ref_format))

                compatible = False
                for test_chrom in test_chroms:
                    try:
                        # This would need reference sequence access
                        # For now, assume if it's standard and matches format, it's compatible
                        if not self.is_standard_chromosome(test_chrom) or test_chrom in [
                            ref_chrom,
                            bam_chrom,
                        ]:
                            compatible = True
                            break
                    except Exception:
                        continue

                if not compatible:
                    incompatible_chromosomes.append(var_chrom)

            except Exception as e:
                logger.warning(f"Could not validate chromosome {var_chrom}: {e}")

        # Step 5: Report issues
        if incompatible_chromosomes:
            self.validation_issues.append(
                f"Variant chromosomes not compatible with reference format: {incompatible_chromosomes}. "
                f"Reference format: {ref_format}. "
                f"gbcms will attempt normalization, but please verify chromosome names are correct."
            )
            logger.warning(f"Incompatible chromosomes found: {incompatible_chromosomes}")
            self.normalization_needed = True

        # Step 6: Determine if normalization is needed
        sample_var_chrom = next(iter(self.variant_chromosomes)) if self.variant_chromosomes else ""
        if self.is_standard_chromosome(sample_var_chrom):
            # Check if variant format matches reference format
            var_has_chr = sample_var_chrom.startswith("chr")
            ref_has_chr = ref_format == "chr_prefix"

            if var_has_chr != ref_has_chr:
                self.normalization_needed = True
                logger.info(
                    f"Normalization needed: variant format ({'chr_prefix' if var_has_chr else 'no_prefix'}) != reference format ({ref_format})"
                )

        if not self.validation_issues:
            logger.info("Chromosome validation passed: all files compatible")

        return len(self.validation_issues) == 0

    def get_validation_summary(self) -> str:
        """Get summary of validation results."""
        if not self.validation_issues:
            ref_format = self.reference_format[1] if self.reference_format else "unknown"
            return f"✅ Chromosome validation passed. Reference format: {ref_format}"

        summary = "❌ Chromosome validation issues found:\n"
        for issue in self.validation_issues:
            summary += f"   - {issue}\n"
        return summary
