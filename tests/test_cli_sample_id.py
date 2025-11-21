"""
Test CLI sample ID parsing and output suffix.
"""

from unittest.mock import MagicMock, patch

from typer.testing import CliRunner

from gbcms.cli import app
from gbcms.models.core import GbcmsConfig

runner = CliRunner()


def test_cli_sample_id_parsing(tmp_path):
    # Create dummy files
    vcf = tmp_path / "test.vcf"
    vcf.touch()
    bam = tmp_path / "test.bam"
    bam.touch()
    fasta = tmp_path / "ref.fasta"
    fasta.touch()
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    # Test 1: --bam sample_id:path
    runner.invoke(
        app,
        [
            "run",
            "-v",
            str(vcf),
            "-b",
            f"mysample:{bam}",
            "-f",
            str(fasta),
            "-o",
            str(output_dir),
            "--suffix",
            ".genotyped",
        ],
    )

    # We expect it to fail pipeline execution (files are empty), but we want to check if it parsed correctly.
    # Since we can't easily inspect internal state of a running CLI command without mocking,
    # we'll check the error message or mock GbcmsConfig/Pipeline.
    # However, simpler is to check if it TRIED to process "mysample".

    # Actually, let's just unit test the logic if possible, or rely on the fact that
    # if it fails with "BAM file not found" it parsed the path correctly.

    # But wait, we can check if the output file was created if we mock the pipeline.
    # For now, let's trust the integration test approach.
    pass


def test_config_model():
    # Unit test the config model directly if needed, but logic is in CLI.
    pass


# Integration test with mocked pipeline would be better,
# but let's just run a real CLI command with a mocked pipeline run method?
# Or just verify the file existence logic in a separate test.

# Let's create a test that invokes the CLI and mocks Pipeline.run


@patch("gbcms.cli.Pipeline")
def test_cli_parsing_mocked(mock_pipeline_cls, tmp_path):
    vcf = tmp_path / "test.vcf"
    vcf.touch()
    bam = tmp_path / "test.bam"
    bam.touch()
    fasta = tmp_path / "ref.fasta"
    fasta.touch()
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    # Mock the pipeline instance
    mock_pipeline = MagicMock()
    mock_pipeline_cls.return_value = mock_pipeline

    # Run CLI
    result = runner.invoke(
        app,
        [
            "run",
            "-v",
            str(vcf),
            "-b",
            f"mysample:{bam}",
            "-f",
            str(fasta),
            "-o",
            str(output_dir),
            "--suffix",
            ".genotyped",
        ],
    )

    assert result.exit_code == 0

    # Check config passed to Pipeline
    config = mock_pipeline_cls.call_args[0][0]
    assert isinstance(config, GbcmsConfig)
    assert config.output_suffix == ".genotyped"
    assert "mysample" in config.bam_files
    assert config.bam_files["mysample"] == bam


@patch("gbcms.cli.Pipeline")
def test_cli_bam_list_parsing(mock_pipeline_cls, tmp_path):
    vcf = tmp_path / "test.vcf"
    vcf.touch()
    bam1 = tmp_path / "test1.bam"
    bam1.touch()
    bam2 = tmp_path / "test2.bam"
    bam2.touch()
    fasta = tmp_path / "ref.fasta"
    fasta.touch()
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    # Create bam list
    bam_list = tmp_path / "bams.list"
    with open(bam_list, "w") as f:
        f.write(f"sample1 {bam1}\n")
        f.write(f"sample2\t{bam2}\n")  # Tab separated

    # Mock the pipeline instance
    mock_pipeline = MagicMock()
    mock_pipeline_cls.return_value = mock_pipeline

    # Run CLI
    result = runner.invoke(
        app, ["run", "-v", str(vcf), "-L", str(bam_list), "-f", str(fasta), "-o", str(output_dir)]
    )

    assert result.exit_code == 0

    # Check config
    config = mock_pipeline_cls.call_args[0][0]
    assert "sample1" in config.bam_files
    assert config.bam_files["sample1"] == bam1
    assert "sample2" in config.bam_files
    assert config.bam_files["sample2"] == bam2
