"""Tests for CLI module."""

import pytest
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()


def test_cli_version():
    """Test --version flag."""
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert "py-gbcms" in result.stdout
    assert "version" in result.stdout


def test_cli_help():
    """Test --help flag."""
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "Python implementation of gbcms" in result.stdout
    assert "--fasta" in result.stdout
    assert "--bam" in result.stdout


def test_cli_missing_required_args():
    """Test CLI with missing required arguments."""
    result = runner.invoke(app, [])
    assert result.exit_code != 0


def test_cli_missing_fasta(sample_bam, sample_vcf, temp_dir):
    """Test CLI with missing FASTA file."""
    result = runner.invoke(
        app,
        [
            "--fasta", str(temp_dir / "nonexistent.fa"),
            "--bam", f"sample1:{sample_bam}",
            "--vcf", str(sample_vcf),
            "--output", str(temp_dir / "output.txt"),
        ],
    )
    assert result.exit_code != 0


def test_cli_mutually_exclusive_maf_vcf(sample_fasta, sample_bam, sample_vcf, sample_maf, temp_dir):
    """Test CLI with both MAF and VCF specified."""
    result = runner.invoke(
        app,
        [
            "--fasta", str(sample_fasta),
            "--bam", f"sample1:{sample_bam}",
            "--vcf", str(sample_vcf),
            "--maf", str(sample_maf),
            "--output", str(temp_dir / "output.txt"),
        ],
    )
    assert result.exit_code != 0
    assert "mutually exclusive" in result.stdout.lower()


def test_cli_no_bam_files(sample_fasta, sample_vcf, temp_dir):
    """Test CLI without BAM files."""
    result = runner.invoke(
        app,
        [
            "--fasta", str(sample_fasta),
            "--vcf", str(sample_vcf),
            "--output", str(temp_dir / "output.txt"),
        ],
    )
    assert result.exit_code != 0


def test_cli_no_variant_files(sample_fasta, sample_bam, temp_dir):
    """Test CLI without variant files."""
    result = runner.invoke(
        app,
        [
            "--fasta", str(sample_fasta),
            "--bam", f"sample1:{sample_bam}",
            "--output", str(temp_dir / "output.txt"),
        ],
    )
    assert result.exit_code != 0


def test_cli_invalid_bam_format(sample_fasta, sample_vcf, temp_dir):
    """Test CLI with invalid BAM format."""
    result = runner.invoke(
        app,
        [
            "--fasta", str(sample_fasta),
            "--bam", "invalid_format",
            "--vcf", str(sample_vcf),
            "--output", str(temp_dir / "output.txt"),
        ],
    )
    assert result.exit_code != 0


def test_cli_valid_vcf_run(sample_fasta, sample_bam, sample_vcf, temp_dir):
    """Test valid CLI run with VCF."""
    output_file = temp_dir / "output.txt"
    result = runner.invoke(
        app,
        [
            "--fasta", str(sample_fasta),
            "--bam", f"sample1:{sample_bam}",
            "--vcf", str(sample_vcf),
            "--output", str(output_file),
        ],
    )
    
    # Should complete successfully
    if result.exit_code != 0:
        print(result.stdout)
        print(result.exception)
    
    assert result.exit_code == 0
    assert output_file.exists()


def test_cli_with_threads(sample_fasta, sample_bam, sample_vcf, temp_dir):
    """Test CLI with multiple threads."""
    output_file = temp_dir / "output.txt"
    result = runner.invoke(
        app,
        [
            "--fasta", str(sample_fasta),
            "--bam", f"sample1:{sample_bam}",
            "--vcf", str(sample_vcf),
            "--output", str(output_file),
            "--thread", "2",
        ],
    )
    
    assert result.exit_code == 0


def test_cli_with_quality_filters(sample_fasta, sample_bam, sample_vcf, temp_dir):
    """Test CLI with quality filters."""
    output_file = temp_dir / "output.txt"
    result = runner.invoke(
        app,
        [
            "--fasta", str(sample_fasta),
            "--bam", f"sample1:{sample_bam}",
            "--vcf", str(sample_vcf),
            "--output", str(output_file),
            "--maq", "30",
            "--baq", "20",
            "--filter-duplicate",
        ],
    )
    
    assert result.exit_code == 0


def test_cli_verbose_mode(sample_fasta, sample_bam, sample_vcf, temp_dir):
    """Test CLI with verbose mode."""
    output_file = temp_dir / "output.txt"
    result = runner.invoke(
        app,
        [
            "--fasta", str(sample_fasta),
            "--bam", f"sample1:{sample_bam}",
            "--vcf", str(sample_vcf),
            "--output", str(output_file),
            "--verbose",
        ],
    )
    
    assert result.exit_code == 0
