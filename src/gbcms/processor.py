"""Main processing logic for GetBaseCounts."""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

import pysam
from rich.progress import BarColumn, Progress, SpinnerColumn, TaskProgressColumn, TextColumn

from .config import Config
from .counter import BaseCounter
from .output import OutputFormatter
from .reference import ReferenceSequence
from .variant import VariantEntry, VariantLoader

logger = logging.getLogger(__name__)


class VariantProcessor:
    """Main processor for counting bases in variants."""

    def __init__(self, config: Config):
        """
        Initialize variant processor.

        Args:
            config: Configuration object
        """
        self.config = config
        self.reference = ReferenceSequence(config.fasta_file)
        self.counter = BaseCounter(config)
        self.sample_order = list(config.bam_files.keys())

    def process(self) -> None:
        """Main processing pipeline."""
        # Load variants
        loader = VariantLoader(reference_getter=self.reference.get_base)
        variants = self._load_all_variants(loader)

        if not variants:
            logger.warning("No variants to process")
            return

        # Sort and index variants
        variants = self._sort_and_index_variants(variants)

        # Initialize counts for all samples
        for variant in variants:
            variant.initialize_counts(self.sample_order)

        # Create variant blocks for parallel processing
        variant_blocks = self._create_variant_blocks(variants)

        logger.info(f"Created {len(variant_blocks)} variant blocks for processing")

        # Process each BAM file
        for sample_name, bam_path in self.config.bam_files.items():
            self._process_bam_file(sample_name, bam_path, variants, variant_blocks)

        # Write output
        self._write_output(variants)

        # Cleanup
        self.reference.close()
        logger.info("Finished processing")

    def _load_all_variants(self, loader: VariantLoader) -> list[VariantEntry]:
        """Load all variants from input files."""
        all_variants = []

        for variant_file in self.config.variant_files:
            if self.config.input_is_maf:
                variants = loader.load_maf(variant_file)
            else:
                variants = loader.load_vcf(variant_file)
            all_variants.extend(variants)

        logger.info(f"Total variants loaded: {len(all_variants)}")
        return all_variants

    def _sort_and_index_variants(self, variants: list[VariantEntry]) -> list[VariantEntry]:
        """Sort variants and identify duplicates."""
        logger.info("Sorting variants")
        variants.sort()

        logger.info("Indexing variants")
        duplicate_map: dict[tuple, VariantEntry] = {}

        for variant in variants:
            key = variant.get_variant_key()
            if key not in duplicate_map:
                duplicate_map[key] = variant
            else:
                # Mark as duplicate
                variant.duplicate_variant_ptr = duplicate_map[key]

        return variants

    def _create_variant_blocks(self, variants: list[VariantEntry]) -> list[tuple[int, int]]:
        """
        Create blocks of variants for parallel processing.

        Returns:
            List of (start_index, end_index) tuples
        """
        if not variants:
            return []

        blocks = []
        start_idx = 0
        current_count = 0

        for i in range(len(variants)):
            current_count += 1

            # Check if we should create a new block
            should_break = False

            if current_count >= self.config.max_block_size:
                should_break = True
            elif i > start_idx:
                # Check chromosome change or distance
                if variants[i].chrom != variants[start_idx].chrom:
                    should_break = True
                elif variants[i].pos - variants[start_idx].pos > self.config.max_block_dist:
                    should_break = True

            if should_break:
                blocks.append((start_idx, i - 1))
                start_idx = i
                current_count = 1

        # Add final block
        blocks.append((start_idx, len(variants) - 1))

        return blocks

    def _process_bam_file(
        self,
        sample_name: str,
        bam_path: str,
        variants: list[VariantEntry],
        variant_blocks: list[tuple[int, int]],
    ) -> None:
        """
        Process a single BAM file.

        Args:
            sample_name: Sample name
            bam_path: Path to BAM file
            variants: List of all variants
            variant_blocks: List of variant block ranges
        """
        logger.info(f"Processing BAM file: {bam_path}")

        if self.config.num_threads > 1:
            self._process_bam_parallel(sample_name, bam_path, variants, variant_blocks)
        else:
            self._process_bam_sequential(sample_name, bam_path, variants, variant_blocks)

    def _process_bam_sequential(
        self,
        sample_name: str,
        bam_path: str,
        variants: list[VariantEntry],
        variant_blocks: list[tuple[int, int]],
    ) -> None:
        """Process BAM file sequentially."""
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
            ) as progress:
                task = progress.add_task(
                    f"[cyan]Processing {sample_name}...", total=len(variant_blocks)
                )

                for start_idx, end_idx in variant_blocks:
                    self._process_variant_block(bam, sample_name, variants, start_idx, end_idx)
                    progress.update(task, advance=1)

    def _process_bam_parallel(
        self,
        sample_name: str,
        bam_path: str,
        variants: list[VariantEntry],
        variant_blocks: list[tuple[int, int]],
    ) -> None:
        """Process BAM file in parallel."""
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
        ) as progress:
            task = progress.add_task(
                f"[cyan]Processing {sample_name}...", total=len(variant_blocks)
            )

            with ThreadPoolExecutor(max_workers=self.config.num_threads) as executor:
                futures = []

                for start_idx, end_idx in variant_blocks:
                    future = executor.submit(
                        self._process_variant_block_thread_safe,
                        bam_path,
                        sample_name,
                        variants,
                        start_idx,
                        end_idx,
                    )
                    futures.append(future)

                for future in as_completed(futures):
                    future.result()  # Raise any exceptions
                    progress.update(task, advance=1)

    def _process_variant_block_thread_safe(
        self,
        bam_path: str,
        sample_name: str,
        variants: list[VariantEntry],
        start_idx: int,
        end_idx: int,
    ) -> None:
        """Process a variant block in a thread-safe manner."""
        # Each thread opens its own BAM file handle
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            self._process_variant_block(bam, sample_name, variants, start_idx, end_idx)

    def _process_variant_block(
        self,
        bam: pysam.AlignmentFile,
        sample_name: str,
        variants: list[VariantEntry],
        start_idx: int,
        end_idx: int,
    ) -> None:
        """
        Process a block of variants.

        Args:
            bam: Open BAM file handle
            sample_name: Sample name
            variants: List of all variants
            start_idx: Start index in variants list
            end_idx: End index in variants list
        """
        start_variant = variants[start_idx]
        end_variant = variants[end_idx]

        # Fetch alignments for the region
        try:
            alignments = list(
                bam.fetch(
                    start_variant.chrom,
                    start_variant.pos,
                    end_variant.pos + 2,  # Buffer for indels
                )
            )
        except Exception as e:
            logger.error(
                f"Error fetching alignments for region "
                f"{start_variant.chrom}:{start_variant.pos}-{end_variant.pos}: {e}"
            )
            return

        # Filter alignments
        filtered_alignments = [aln for aln in alignments if not self.counter.filter_alignment(aln)]

        # Process each variant in the block
        for i in range(start_idx, end_idx + 1):
            variant = variants[i]

            # Skip if this is a duplicate variant
            if variant.duplicate_variant_ptr is not None:
                continue

            # Count bases for this variant
            self.counter.count_variant(variant, filtered_alignments, sample_name)

    def _write_output(self, variants: list[VariantEntry]) -> None:
        """Write output file."""
        formatter = OutputFormatter(self.config, self.sample_order)

        if self.config.input_is_maf:
            if self.config.output_maf:
                formatter.write_maf_output(variants)
            else:
                formatter.write_fillout_output(variants)
        else:
            formatter.write_vcf_output(variants)
